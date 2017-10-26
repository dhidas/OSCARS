echo this is include...
echo %INCLUDE%

echo Downloading CUDA toolkit 8
appveyor DownloadFile  https://developer.nvidia.com/compute/cuda/8.0/prod/local_installers/cuda_8.0.44_windows-exe -FileName cuda_8.0.44_windows.exe
echo Installing CUDA toolkit 8
cuda_8.0.44_windows.exe -s compiler_8.0 ^
                           cublas_8.0 ^
                           cublas_dev_8.0 ^
                           cudart_8.0 ^
                           curand_8.0 ^
                           curand_dev_8.0

if NOT EXIST "%ProgramFiles%\NVIDIA GPU Computing Toolkit\CUDA\v8.0\bin\cudart64_80.dll" ( 
echo "Failed to install CUDA"
exit /B 1
)

echo %PATH%


set PATH=%ProgramFiles%\NVIDIA GPU Computing Toolkit\CUDA\v8.0\bin;%ProgramFiles%\NVIDIA GPU Computing Toolkit\CUDA\v8.0\libnvvp;%PATH%
set PATH=C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin;%PATH%
set INCLUDE=C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\INCLUDE;%INCLUDE%
set LIB=C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\LIB;%LIB%
set INCLUDE=C:\Program Files (x86)\Windows Kits\10\Include\10.0.10586\ucrt;%INCLUDE%
set INC=C:\Program Files (x86)\Windows Kits\10\Include\10.0.10586\ucrt;%INCLUDE%
set LIB=C:\Program Files (x86)\Windows Kits\10\Lib\10.0.10586\um\x64;C;\Program Files (x86)\Windows Kits\10\Lib\10.0.10586\ucrt\x64;%LIB%

echo include dir
dir "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\INCLUDE"

nvcc -V

ECHO nvcc -DCUDA -cudart static -shared -Iinclude -I%INC% -c src\OSCARSSR_Cuda.cu -o lib\OSCARSSR_Cuda.o
ECHO nvcc -DCUDA -cudart static -shared -Iinclude -I%INC% -c src\OSCARSTH_Cuda.cu -o lib\OSCARSTH_Cuda.o

dir "lib"
dir "."

echo done build cuda
