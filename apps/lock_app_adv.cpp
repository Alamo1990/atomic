#include "../include/locked_buffer.h"
#include <thread>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <utility>
#include <fstream>
#include <chrono>

#define NTHREADS 4
constexpr int ImageHeight = 36;
constexpr int ImageWidth = 36;
using namespace std;
mutex mtx;

pair<int,vector<vector<unsigned char> > > mandelbrot(double MaxRe, double MinRe, double MinIm, int order){
        vector<vector<unsigned char> > image(ImageHeight, vector<unsigned char>(ImageWidth,0));

    //   double MinRe = -2.0;
    //   double MaxRe = p;
    //   double MinIm = -1.2;
    double MaxIm = MinIm+(MaxRe-MinRe)*ImageHeight/ImageWidth;
    double Re_factor = (MaxRe-MinRe)/(ImageWidth-1);
    double Im_factor = (MaxIm-MinIm)/(ImageHeight-1);
    unsigned MaxIterations = 30;

    for(unsigned y=0; y<ImageHeight; ++y)
    {
        double c_im = MaxIm - y*Im_factor;
        for(unsigned x=0; x<ImageWidth; ++x)
        {
            double c_re = MinRe + x*Re_factor;
            double Z_re = c_re, Z_im = c_im;
            bool isInside = true;
            for(unsigned n=0; n<MaxIterations; ++n)
            {
                double Z_re2 = Z_re*Z_re, Z_im2 = Z_im*Z_im;
                if(Z_re2 + Z_im2 > 4)
                {
                    isInside = false;
                    break;
                }
                Z_im = 2*Z_re*Z_im + c_im;
                Z_re = Z_re2 - Z_im2 + c_re;
            }
            if(isInside) { image[x][y] = 255; }
        }
    }
    return make_pair(order,image);
}

pair<int,vector<vector<complex<double> > > > FFT( pair<int,vector<vector<unsigned char> > > & image ){
    vector<vector<complex<double> > > F(ImageHeight, vector<complex<double> >(ImageWidth));
    complex<double> J (0,1);
    for(int k = 0; k< ImageHeight; k++) {
        for(int l=0; l<ImageWidth; ++l) {
            complex<double> sum = 0;
            for(int m = 0; m< ImageHeight; m++) {
                for(int n=0; n<ImageWidth; ++n) {
                    sum += (image.second[m][n]*1.0) * exp(-1.0*J*2.0*M_PI*((k*m*1.0)/(ImageHeight*1.0) + (l*n)/(ImageWidth*1.0)   ));
                }
            }
            F[k][l] = 1.0/(ImageHeight*ImageWidth) * sum;
        }
    }
    return make_pair(image.first,F);
}


pair<int,vector<vector<unsigned char> > > IFFT(pair<int, vector<vector<complex<double> > > > & image){
    vector<vector<unsigned char > > f(ImageHeight, vector<unsigned char>(ImageWidth));
    complex<double> J (0,1);
    for(int m = 0; m< ImageHeight; ++m) {
        for(int n=0; n<ImageWidth; ++n) {
            complex<double> sum = 0;
            for(int k = 0; k< ImageHeight; ++k) {
                for(int l=0; l<ImageWidth; ++l) {
                        sum += (image.second[k][l]) * exp(J*2.0*M_PI*((k*m)/(ImageHeight*1.0) + (l*n)/(ImageWidth*1.0) ));
                }
            }
            f[m][n] = (unsigned char)round(sum.real());
        }
    }
    return make_pair(image.first, f);
}

pair<int, vector< vector< complex<double> > > > Blur( pair<int, vector<vector<complex<double> > > > & image){
    vector<vector<complex<double> > > F(ImageHeight, vector<complex<double> >(ImageWidth, complex<double>(0,0)));
    constexpr double kernel[3][3] = {{1/18.0, 1/18.0, 1/18.0},{1/18.0, 10/18.0, 1/18.0}, {1/18.0, 1/18.0, 1/18.0} };
    for(int y=1; y<ImageHeight-1; y++) {
        for(int x=1; x<ImageWidth-1; x++) {
            complex<double> sum (0,0);
            for(int k=-1; k<= 1; k++) {
                for(int j = -1; j<=1; j++) {
                        sum += kernel[j+1][k+1] * image.second[y-j][x-k];
                }
            }
            F[y][x] = sum;
        }
    }
    return make_pair(image.first,F);
}


void print(pair<int, vector<vector<unsigned char> > > image){
    ofstream file;
    file.open("output/"+to_string(image.first)+".txt");
    for(int i = 0; i<ImageHeight; i++) {
        for(int j=0; j<ImageWidth; j++) {
                file<<(int)image.second[i][j];
        }
        file<<"\n";
    }
    file.close();
}

using uchar_matrix      =   vector<vector<unsigned char>>;
using complex_matrix    =   vector<vector<complex<double>>>;

void scheduler(int nthreads,
    int nitems,
    locked_buffer<pair<int, uchar_matrix > >*  queue0,
    vector< unique_ptr< locked_buffer<pair<int,uchar_matrix > > > > & queues1){

    vector<int> sent(nthreads);
    vector<int> tosend(nthreads, nitems/nthreads-1);
    for (int turn=0; turn<nthreads; turn++)
        if(nitems%nthreads>turn)
            tosend.at(turn)++;
    for (int i = 0; i < nitems; i++) {
            int turn = i%nthreads;

            auto image = get<1>(queue0->get());
            queues1[turn]->put(image, tosend.at(turn)==sent.at(turn));
            sent.at(turn)++;

    }
}

void looper(int mode,
    int nitems,
    locked_buffer<pair<int, uchar_matrix > >*                   queue0,
    unique_ptr<locked_buffer<pair<int, uchar_matrix > > >&      queue1,
    unique_ptr< locked_buffer<pair<int, complex_matrix > > >&   queue2,
    unique_ptr< locked_buffer<pair<int, complex_matrix > > >&   queue3,
    locked_buffer<pair<int, uchar_matrix > >*                   queue4
            ){

        enum exec { Mandelbrot_fn, FFT_fn, Blur_fn, IFFT_fn, Print_fn };
        bool last = false;
        double MaxRe;
        double MinRe;
        double MinIm;

        switch (mode) {
        case Mandelbrot_fn:
            for(int i = 0; i<nitems; i++)
            {
                MaxRe = 0.1 + i *0.1;
                MinRe = -2.0 - i *0.1;
                MinIm = -1.2 - i *0.1;
                queue0->put(mandelbrot(MaxRe, MinRe, MinIm, i), i == (nitems-1));
            } break;
        case FFT_fn:
            while(!last)
            {
                auto elem = queue1->get();
                last = get<0>(elem);
                auto image = get<1>(elem);
                queue2->put(FFT(image), last);
            } break;
        case Blur_fn:
            while(!last)
            {
                auto elem = queue2->get();
                last = get<0>(elem);
                auto image = get<1>(elem);
                queue3->put(Blur(image), last);
            } break;
        case IFFT_fn:
            while(!last)
            {
                auto elem = queue3->get();
                last = get<0>(elem);
                auto image = get<1>(elem);
                queue4->put(IFFT(image), last);
            } break;
        case Print_fn:
            for(int i = 0; i<nitems; i++)
            {
                auto image = get<1>(queue4->get());
                print(image);
            } break;
        }
}

int main(int argc, char* argv[]){
        chrono::time_point<chrono::system_clock> start, end;
        start = chrono::system_clock::now();
        if(argc != 4) {
                cerr<<"Wrong arguments"<<endl;
                cerr<<"Valid formats: "<<endl;
                cerr<< " " << argv[0] << "nitems" << endl;
        }
        const long nitems = stol(argv[1]);
        const int buff_size = stoi(argv[2]);
        const int nthreads = stoi(argv[3]);


        auto queue0 = new locked_buffer<pair<int,vector<vector<unsigned char> > > >(buff_size);
        vector< unique_ptr< locked_buffer<pair<int, uchar_matrix > > > >    queues1;
        vector< unique_ptr< locked_buffer<pair<int, complex_matrix > > > >  queues2;
        vector< unique_ptr< locked_buffer<pair<int, complex_matrix > > > >  queues3;
        for(int i=0; i<nthreads; i++) {
                queues1.push_back( unique_ptr< locked_buffer<pair<int, uchar_matrix > > >   (new locked_buffer<pair<int, uchar_matrix > >(buff_size)));
                queues2.push_back( unique_ptr< locked_buffer<pair<int, complex_matrix > > > (new locked_buffer<pair<int, complex_matrix > >(buff_size)));
                queues3.push_back( unique_ptr< locked_buffer<pair<int, complex_matrix > > > (new locked_buffer<pair<int, complex_matrix > >(buff_size)));
        }
        locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue4 = new locked_buffer<pair<int,vector<vector<unsigned char> > > >(buff_size);

        vector<thread> threads(3+nthreads*3);

        threads.at(0) = thread(looper, 0, nitems, queue0, ref(queues1[0]), ref(queues2[0]), ref(queues3[0]), queue4);
        threads.at(1) = thread(scheduler, nthreads, nitems, queue0, ref(queues1));
        for (int i = 0; i < nthreads; i++) {
                threads.at(3+i*3) = thread(looper, 1, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                threads.at(4+i*3) = thread(looper, 2, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                threads.at(5+i*3) = thread(looper, 3, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                threads.at(3+i*3).detach();
                threads.at(4+i*3).detach();
                threads.at(5+i*3).detach();
        }
        threads.at(2+nthreads*3) = thread(looper, 4, nitems, queue0, ref(queues1[0]), ref(queues2[0]), ref(queues3[0]), queue4);


        threads.at(0).join();
        threads.at(1).join();
        threads.at(2+nthreads*3).join();

        end = chrono::system_clock::now();

        int elapsed_milliseconds = chrono::duration_cast<chrono::milliseconds>
                                           (end-start).count();
        time_t end_time = chrono::system_clock::to_time_t(end);
        cout << "End time:  " << ctime(&end_time)
                  << "\ntime elapsed: " << elapsed_milliseconds << "ms\n";

        return 0;
}
