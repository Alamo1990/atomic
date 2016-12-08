#include "../include/atomic_buffer.h"
#include <thread>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <utility>
#include <fstream>
#include <atomic>

#define NTHREADS 4
using namespace std;
constexpr int ImageHeight = 36;
constexpr int ImageWidth = 36;
atomic_flag atomFlag = ATOMIC_FLAG_INIT;

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

void scheduler(int nthreads, int nitems, atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >* queue0,
               std::vector< std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > > > >&  queues1){
        for (int i = 0; i < nitems; i++) {
                int turn = i%nthreads;
                auto image = std::get<1>(queue0->get());
                queues1[turn]->put(image, i==(nitems-1));
        }
}

void looper(int mode,
            int nitems,
            atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >* queue0,
            std::unique_ptr<atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > > >& queue1,
            std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<std::complex<double> > > > > >& queue2,
            std::unique_ptr< atomic_buffer<std::pair<int, std::vector< std::vector< std::complex<double> > > > > >& queue3,
            atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >* queue4
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
                        queue0->put(mandelbrot(MaxRe, MinRe, MinIm, i), last); // no se si true o false
                        break;
                case F: {
                        auto image = std::get<1>(queue1->get());
                        queue2->put(FFT(image), last);
                }
                break;
                case B: {
                        auto image = std::get<1>(queue2->get());
                        queue3->put(Blur(image), last);
                        break;
                }
                case I: {
                        auto image = std::get<1>(queue3->get());
                        queue4->put(IFFT(image), last);
                } break;
                case P:
                        while(atomFlag.test_and_set());
                        print(std::get<1>(queue4->get()));
                        atomFlag.clear();
                        break;
                }
        }
}

int main(int argc, char* argv[]){
        if(argc != 4) {
                std::cerr<<"Wrong arguments"<<std::endl;
                std::cerr<<"Valid formats: "<<std::endl;
                std::cerr<< " " << argv[0] << "nitems" << std::endl;
        }
        const long nitems = std::stol(argv[1]);
        const int buff_size = std::stoi(argv[2]);
        const int nthreads = std::stoi(argv[3]);


        atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >* queue0 = new atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >(buff_size);

        std::vector< std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > > > >  queues1;
        std::vector< std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<std::complex<double> > > > > > >  queues2;
        std::vector< std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<std::complex<double> > > > > > >  queues3;
        for(int i=0; i<nthreads; i++) {
                queues1.push_back( std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > > > (new atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >(buff_size)));
                queues2.push_back( std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<std::complex<double> > > > > > (new atomic_buffer<std::pair<int,std::vector<std::vector<std::complex<double> > > > >(buff_size)));
                queues3.push_back( std::unique_ptr< atomic_buffer<std::pair<int,std::vector<std::vector<std::complex<double> > > > > > (new atomic_buffer<std::pair<int,std::vector<std::vector<std::complex<double> > > > >(buff_size)));
        }
        atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >* queue4 = new atomic_buffer<std::pair<int,std::vector<std::vector<unsigned char> > > >(buff_size);


        std::vector<std::thread> threads(2+nthreads*3+1);

        threads.at(0) = std::thread(looper, 0, nitems, queue0, ref(queues1[0]), ref(queues2[0]), ref(queues3[0]), queue4);
        std::cout << "aquiii1" << std::endl;
        threads.at(1) = std::thread(scheduler, nthreads, nitems, queue0, ref(queues1));
        std::cout << "aquiii2" << std::endl;
        for (int i = 0; i < nthreads; i++) {
                threads.at(3+i*3) = std::thread(looper, 1, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                std::cout << "aquiii3" << std::endl;
                threads.at(4+i*3) = std::thread(looper, 2, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                std::cout << "aquiii4" << std::endl;
                threads.at(5+i*3) = std::thread(looper, 3, nitems, queue0, ref(queues1[i]), ref(queues2[i]), ref(queues3[i]), queue4);
                std::cout << "aquiii5" << std::endl;
        }
        threads.at(2+nthreads*3) = std::thread(looper, 4,nitems, queue0, ref(queues1[0]), ref(queues2[0]), ref(queues3[0]), queue4);

        std::cout << "aquiii6" << std::endl;

        for(auto& th : threads) th.join();

        return 0;
}
