#include "../include/locked_buffer.h"
#include <thread>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <utility>
#include <fstream>

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

void scheduler(int nthreads, int nitems, locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue0,
               vector< unique_ptr< locked_buffer<pair<int,vector<vector<unsigned char> > > > > > & queues1){
        cout << "schestart" << endl;
        for (int i = 0; i < nitems; i++) {
                cout << "schloopin: "<<i << endl;
                int turn = i%nthreads;
                cout << "schget: "<< i << endl;
                auto image = get<1>(queue0->get());
                cout << "schpreput: "<< i << endl;
                queues1[turn]->put(image, i==(nitems-1));
        }
}

void looper(int mode,
            int nitems,
            locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue0,
            unique_ptr<locked_buffer<pair<int,vector<vector<unsigned char> > > > >& queue1,
            unique_ptr< locked_buffer<pair<int,vector<vector<complex<double> > > > > >& queue2,
            unique_ptr< locked_buffer<pair<int, vector< vector< complex<double> > > > > >& queue3,
            locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue4
            ){
        enum exec { M, F, B, I, P };
        bool last = false;
        for(int i = 0; i<nitems; i++) {
                double MaxRe = 0.1 + i *0.1;
                double MinRe = -2.0 - i *0.1;
                double MinIm = -1.2 - i *0.1;
                if(i == (nitems-1)) last=true;

                switch (mode) {
                case M:
                        cout << "prebrot" << endl;
                        queue0->put(mandelbrot(MaxRe, MinRe, MinIm, i), last); // no se si true o false
                        cout << "postbrot" << endl;
                        break;
                case F: {
                        cout << "getFFT" << endl;
                        auto image = get<1>(queue1->get());
                        cout << "preFFT" << endl;
                        queue2->put(FFT(image), last);
                        cout << "postFFT" << endl;
                }
                break;
                case B: {
                        auto image = get<1>(queue2->get());
                        cout << "preblur" << endl;
                        queue3->put(Blur(image), last);
                        cout << "postblur" << endl;
                        break;
                }
                case I: {
                        auto image = get<1>(queue3->get());
                        cout << "preinv" << endl;
                        queue4->put(IFFT(image), last);
                        cout << "postinv" << endl;
                } break;
                case P: {
                        mtx.lock();
                        auto image = get<1>(queue4->get());
                        mtx.unlock();
                        print(image);
                        break;
                }
                }
        }
}
int main(int argc, char* argv[]){
        if(argc != 4) {
                cerr<<"Wrong arguments"<<endl;
                cerr<<"Valid formats: "<<endl;
                cerr<< " " << argv[0] << "nitems" << endl;
        }
        const long nitems = stol(argv[1]);
        const int buff_size = stoi(argv[2]);
        const int nthreads = stoi(argv[3]);


        locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue0 = new locked_buffer<pair<int,vector<vector<unsigned char> > > >(buff_size);

        vector< unique_ptr< locked_buffer<pair<int,vector<vector<unsigned char> > > > > >  queues1;
        vector< unique_ptr< locked_buffer<pair<int,vector<vector<complex<double> > > > > > >  queues2;
        vector< unique_ptr< locked_buffer<pair<int,vector<vector<complex<double> > > > > > >  queues3;
        for(int i=0; i<nthreads; i++) {
                queues1.push_back( unique_ptr< locked_buffer<pair<int,vector<vector<unsigned char> > > > > (new locked_buffer<pair<int,vector<vector<unsigned char> > > >(buff_size)));
                queues2.push_back( unique_ptr< locked_buffer<pair<int,vector<vector<complex<double> > > > > > (new locked_buffer<pair<int,vector<vector<complex<double> > > > >(buff_size)));
                queues3.push_back( unique_ptr< locked_buffer<pair<int,vector<vector<complex<double> > > > > > (new locked_buffer<pair<int,vector<vector<complex<double> > > > >(buff_size)));
        }
        locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue4 = new locked_buffer<pair<int,vector<vector<unsigned char> > > >(buff_size);


        vector<thread> threads(2+nthreads*3+1);

        threads.at(0) = thread(looper, 0, nitems, queue0, ref(queues1[0]), ref(queues2[0]), ref(queues3[0]), queue4);
        cout << "0" << endl;
        threads.at(1) = thread(scheduler, nthreads, nitems, queue0, ref(queues1));
        cout << "1" << endl;
        for (int i = 0; i < nthreads; i++) {
                cout << "loopstart: "<< i << endl;
                threads.at(3+i*3) = thread(looper, 1, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                cout << "2: "<< i << endl;
                threads.at(4+i*3) = thread(looper, 2, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                cout << "3: "<<i << endl;
                threads.at(5+i*3) = thread(looper, 3, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                cout << "4:"<< i << endl;
        }
        threads.at(2+nthreads*3) = thread(looper, 4,nitems, queue0, ref(queues1[0]), ref(queues2[0]), ref(queues3[0]), queue4);
        cout << "5" << endl;

        for(auto& th : threads) th.join();

        return 0;
}
