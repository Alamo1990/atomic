#include "../include/locked_buffer.h"
#include <thread>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <utility>
#include <fstream>
#include <algorithm>

#define NTHREADS 4
constexpr int ImageHeight = 36;
constexpr int ImageWidth = 36;
using namespace std;
mutex mtx;

std::pair<int,std::vector<std::vector<unsigned char> > > mandelbrot(double MaxRe, double MinRe, double MinIm, int order){
        std::vector<std::vector<unsigned char> > image(ImageHeight, std::vector<unsigned char>(ImageWidth,0));

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

void scheduler(int nthreads, std::vector<int> assignations, int nitems, locked_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >* queue0,
               std::vector< std::unique_ptr< locked_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > > > > & queues1){
        // std::cout << "schestart" << std::endl;
        int threadTurn=0;
        while(std::all_of(assignations.begin(), assignations.end(), [] (int i) { return i==0; })) {
                if(assignations.at(threadTurn)>0) {
                        auto image = get<1>(queue0->get());
                        // cout << "schpreput: "<< i << endl;
                        queues1[threadTurn]->put(image, assignations.at(threadTurn)==(nitems-1));
                        assignations.at(threadTurn)--;
                }
                threadTurn = (threadTurn+1) % nthreads;
        }

        for (int i = 0; i < nitems; i++) {
                // cout << "schloopin: "<<i << endl;
                int turn = i%nthreads;

                // cout << "schget: "<< i << endl;
                auto image = get<1>(queue0->get());
                // cout << "schpreput: "<< i << endl;
                queues1[turn]->put(image, i==(nitems-1));
        }

}

void looper(int mode, std::vector<int> assignations, int thID, int nitems,
            locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue0,
            unique_ptr<locked_buffer<pair<int,vector<vector<unsigned char> > > > >& queue1,
            unique_ptr< locked_buffer<pair<int,vector<vector<complex<double> > > > > >& queue2,
            unique_ptr< locked_buffer<pair<int, vector< vector< complex<double> > > > > >& queue3,
            locked_buffer<pair<int,vector<vector<unsigned char> > > >* queue4
            ){
        enum exec { M, F, B, I, P };
        bool last = false;

        double MaxRe;
        double MinRe;
        double MinIm;

        int assigned = assignations.at(thID);

        std::cout << "inThread -> mode: " << mode << " | thID: " << thID << " | assigned: " << assigned << std::endl;

        switch (mode) {
        case M:
                for(int i = 0; i<nitems; i++) {
                        MaxRe = 0.1 + i *0.1;
                        MinRe = -2.0 - i *0.1;
                        MinIm = -1.2 - i *0.1;
                        if(i == (nitems-1)) last=true;

                        // std::cout << "prebrot: "<< i  << std::endl;
                        queue0->put(mandelbrot(MaxRe, MinRe, MinIm, i), last);
                        // std::cout << "postbrot: "<< i << std::endl;
                } break;
        case F:
                for(int i=0; i<assigned; i++) {
                        if(i == (assigned-1)) last=true;

                        // std::cout << "getFFT: "<< i << std::endl;
                        auto image = std::get<1>(queue1->get());
                        // std::cout << "preFFT: "<< i << std::endl;
                        queue2->put(FFT(image), last);
                        // std::cout << "postFFT: "<< i << std::endl;
                } break;
        case B:
                for(int i=0; i<assigned; i++) {
                        if(i == (assigned-1)) last=true;

                        auto image = std::get<1>(queue2->get());
                        // std::cout << "preblur: "<< i << std::endl;
                        queue3->put(Blur(image), last);
                        // std::cout << "postblur: "<< i << std::endl;
                } break;
        case I:
                for(int i=0; i<assigned; i++) {
                        if(i == (assigned-1)) last=true;

                        auto image = std::get<1>(queue3->get());
                        // std::cout << "preinv: "<< i  << std::endl;
                        queue4->put(IFFT(image), last);
                        // std::cout << "postinv: "<< i  << std::endl;
                } break;
        case P:
                for(int i = 0; i<nitems; i++) {
                        if(i == (nitems-1)) last=true;

                        // std::cout << "printpreget: "<< i  << std::endl;
                        auto image = std::get<1>(queue4->get());
                        // std::cout << "itermediateprint: "<< i  << std::endl;
                        print(image);
                        // std::cout << "postprint: "<< i  << std::endl;
                } break;
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

        std::vector<int> assignations(nthreads, nitems/nthreads);
        for(int i=0; i<nitems%nthreads; i++) {
                assignations.at(i)++;
        }

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

        std::vector<std::thread> threads(3+nthreads*3);

        threads.at(0) = std::thread(looper, 0, assignations, 0, nitems, queue0, std::ref(queues1[0]), std::ref(queues2[0]), std::ref(queues3[0]), queue4);
        // std::cout << "0" << std::endl;
        threads.at(1) = std::thread(scheduler, nthreads, assignations, nitems, queue0, std::ref(queues1));
        // std::cout << "1" << std::endl;
        for (int i = 0; i < nthreads; i++) {
                // std::cout << "loopstart: "<< i << std::endl;
                threads.at(3+i*3) = std::thread(looper, 1, assignations, i, nitems, queue0, std::ref(queues1[i]), std::ref(queues2[i]), std::ref(queues3[i]), queue4);
                // std::cout << "2: "<< i << std::endl;
                threads.at(4+i*3) = std::thread(looper, 2, assignations, i, nitems, queue0, std::ref(queues1[i]), std::ref(queues2[i]), std::ref(queues3[i]), queue4);
                // std::cout << "3: "<<i << std::endl;
                threads.at(5+i*3) = std::thread(looper, 3, assignations, i, nitems, queue0, std::ref(queues1[i]), std::ref(queues2[i]), std::ref(queues3[i]), queue4);
                // std::cout << "4:"<< i << std::endl;
                threads.at(3+i*3).detach();
                threads.at(4+i*3).detach();
                threads.at(5+i*3).detach();
        }
        // std::cout << "4,5" << std::endl;
        threads.at(2+nthreads*3) = std::thread(looper, 4, assignations, 0, nitems, queue0, std::ref(queues1[0]), std::ref(queues2[0]), std::ref(queues3[0]), queue4);
        // std::cout << "5" << std::endl;


        threads.at(0).join();
        threads.at(1).join();
        threads.at(2+nthreads*3).join();

        return 0;
}
