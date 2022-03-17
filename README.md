# SVHSS
The code of SVHSS in the paper
## IDE
In the coding, we chose CLion as IDE. Therefore, if you use CLion as IDE it will help you deploy this project mor conveniently.
## Dependencies
### Install m4
    wget http://mirrors.kernel.org/gnu/m4/m4-1.4.13.tar.gz
    tar -xzvf m4-1.4.13.tar.gz
    cd m4-1.4.13
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
### Install autoconf
    wget http://mirrors.kernel.org/gnu/autoconf/autoconf-2.65.tar.gz
    tar -xzvf autoconf-2.65.tar.gz
    cd autoconf-2.65
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
### Install automake
    wget http://mirrors.kernel.org/gnu/automake/automake-1.11.tar.gz
    tar xzvf automake-1.11.tar.gz
    cd automake-1.11
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
### Install libtool
    wget http://mirrors.kernel.org/gnu/libtool/libtool-2.2.6b.tar.gz
    tar xzvf libtool-2.2.6b.tar.gz
    cd libtool-2.2.6b
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
### Install GMP
    tar -jxvf  gmp-6.1.2.tar.bz2
    cd gmp-6.1.2
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
    make check
### Install MPFR
    tar -xzvf mpfr-4.0.1.tar.gz
    cd mpfr-4.0.1
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
### Install FLINT
    tar -xzvf flint-2.5.2.tar.gz
    cd flint-2.5.2
    ./configure --with-gmp  --with-mpfr --prefix=/usr/local
    sudo make
    sudo make install
### Install CMake
    tar -xzvf cmake-3.10.2.tar.gz
    cd cmake-3.10.2
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
