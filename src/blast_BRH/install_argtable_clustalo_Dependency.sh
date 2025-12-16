## argtable was missing
## on macos brew install argtable
## on hpc, download executeable from here https://sourceforge.net/projects/argtable/ (argtable2-x)
## install instructions https://github.com/jonathanmarvens/argtable2

## use this configure for clustalo to tell it where the local install of argtable is
./configure --prefix=/proj/naiss2023-6-65/Milena/software_install/clustal_omega/clustal-omega-1.2.4 CPPFLAGS="-I/proj/naiss2023-6-65/Milena/software_install/clustal_omega/argtable2-13/include" CFLAGS="-I/proj/naiss2023-6-65/Milena/software_install/clustal_omega/argtable2-13/include" LDFLAGS="-L/proj/naiss2023-6-65/Milena/software_install/clustal_omega/argtable2-13/lib"
make
make install

## chatgpt also says this
export LD_LIBRARY_PATH=/proj/naiss2023-6-65/Milena/software_install/clustal_omega/argtable2-13/lib:$LD_LIBRARY_PATH

## then clustalo is installed and in the local bin