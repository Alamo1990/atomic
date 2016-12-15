#include "../include/atomic_buffer.h"
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
        return std::make_pair(order,image);
}

std::pair<int,std::vector<std::vector<std::complex<double> > > > FFT( std::pair<int,std::vector<std::vector<unsigned char> > > & image ){
        std::vector<std::vector<std::complex<double> > > F(ImageHeight, std::vector<std::complex<double> >(ImageWidth));
        std::complex<double> J (0,1);
        for(int k = 0; k< ImageHeight; k++) {
                for(int l=0; l<ImageWidth; ++l) {
                        std::complex<double> sum = 0;
                        for(int m = 0; m< ImageHeight; m++) {
                                for(int n=0; n<ImageWidth; ++n) {
                                        sum += (image.second[m][n]*1.0) * exp(-1.0*J*2.0*M_PI*((k*m*1.0)/(ImageHeight*1.0) + (l*n)/(ImageWidth*1.0)   ));
                                }
                        }

                        F[k][l] = 1.0/(ImageHeight*ImageWidth) * sum;
                }
        }
        return std::make_pair(image.first,F);

}


std::pair<int,std::vector<std::vector<unsigned char> > > IFFT(std::pair<int, std::vector<std::vector<std::complex<double> > > > & image){
        std::vector<std::vector<unsigned char > > f(ImageHeight, std::vector<unsigned char>(ImageWidth));
        std::complex<double> J (0,1);
        for(int m = 0; m< ImageHeight; ++m) {
                for(int n=0; n<ImageWidth; ++n) {
                        std::complex<double> sum = 0;
                        for(int k = 0; k< ImageHeight; ++k) {
                                for(int l=0; l<ImageWidth; ++l) {
                                        sum += (image.second[k][l]) * exp(J*2.0*M_PI*((k*m)/(ImageHeight*1.0) + (l*n)/(ImageWidth*1.0) ));
                                }
                        }
                        f[m][n] = (unsigned char)round(sum.real());
                }
        }

        return std::make_pair(image.first, f);
}

std::pair<int, std::vector< std::vector< std::complex<double> > > > Blur( std::pair<int, std::vector<std::vector<std::complex<double> > > > & image){
        std::vector<std::vector<std::complex<double> > > F(ImageHeight, std::vector<std::complex<double> >(ImageWidth, std::complex<double>(0,0)));
        constexpr double kernel[3][3] = {{1/18.0, 1/18.0, 1/18.0},{1/18.0, 10/18.0, 1/18.0}, {1/18.0, 1/18.0, 1/18.0} };
        for(int y=1; y<ImageHeight-1; y++) {
                for(int x=1; x<ImageWidth-1; x++) {
                        std::complex<double> sum (0,0);
                        for(int k=-1; k<= 1; k++) {
                                for(int j = -1; j<=1; j++) {
                                        sum += kernel[j+1][k+1] * image.second[y-j][x-k];
                                }
                        }
                        F[y][x] = sum;
                }
        }
        return std::make_pair(image.first,F);
}


void print(std::pair<int, std::vector<std::vector<unsigned char> > > image){
        std::ofstream file;
        file.open("output/"+std::to_string(image.first)+".txt");
        for(int i = 0; i<ImageHeight; i++) {
                for(int j=0; j<ImageWidth; j++) {
                        file<<(int)image.second[i][j];
                }
                file<<"\n";
        }
        file.close();
}


//using uchar_matrix = vector<vector<unsigned char> >;
using uchar_matrix      =   vector<vector<unsigned char>>;
using complex_matrix    =   vector<vector<complex<double>>>;

void looper(int mode,
    int nitems,
    atomic_buffer<pair<int, uchar_matrix > >*   queue1,
    atomic_buffer<pair<int, complex_matrix > >* queue2,
    atomic_buffer<pair<int, complex_matrix > >* queue3,
    atomic_buffer<pair<int, uchar_matrix > >*   queue4 ){

    enum exec { Mandelbrot_fn, FFT_fn, Blur_fn, IFFT_fn, Print_fn };
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
            queue1->put(mandelbrot(MaxRe, MinRe, MinIm, i), false);
        } break;
    case FFT_fn:
        for(int i = 0; i<nitems; i++)
        {
            auto image = std::get<1>(queue1->get());
            queue2->put(FFT(image), false);
        } break;
    case Blur_fn:
        for(int i = 0; i<nitems; i++)
        {
            auto image = std::get<1>(queue2->get());
            queue3->put(Blur(image), false);
        } break;
    case IFFT_fn:
        for(int i = 0; i<nitems; i++)
        {
            auto image = std::get<1>(queue3->get());
            queue4->put(IFFT(image), false);
        } break;
    case Print_fn:
        for(int i = 0; i<nitems; i++)
        {
            print(std::get<1>(queue4->get()));
        } break;
    }
}

int main(int argc, char* argv[]){
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        if(argc != 3) {
                cerr<<"Wrong arguments"<<endl;
                cerr<<"Valid formats: "<<endl;
                cerr<< " " << argv[0] << "nitems" << endl;
        }
        const long nitems = stol(argv[1]);
        const int buff_size = stoi(argv[2]);

        auto queue1 = new atomic_buffer<pair<int,uchar_matrix   > >(buff_size);
        auto queue2 = new atomic_buffer<pair<int,complex_matrix > >(buff_size);
        auto queue3 = new atomic_buffer<pair<int,complex_matrix > >(buff_size);
        auto queue4 = new atomic_buffer<pair<int,uchar_matrix   > >(buff_size);

        thread threads[5];

        threads[0] = thread(looper, 0, nitems, queue1, queue2, queue3, queue4);
        threads[1] = thread(looper, 1, nitems, queue1, queue2, queue3, queue4);
        threads[2] = thread(looper, 2, nitems, queue1, queue2, queue3, queue4);
        threads[3] = thread(looper, 3, nitems, queue1, queue2, queue3, queue4);
        threads[4] = thread(looper, 4, nitems, queue1, queue2, queue3, queue4);

        for(auto& th : threads) th.join();

        end = std::chrono::system_clock::now();

        int elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>
                                           (end-start).count();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::cout << "End time:  " << std::ctime(&end_time)
                  << "\ntime elapsed: " << elapsed_milliseconds << "ms\n";

        return 0;
}
