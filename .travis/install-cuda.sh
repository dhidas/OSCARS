export CUDA_VERSION=7.0-28


wget "http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1204/x86_64/cuda-repo-ubuntu1204_${CUDA_VERSION}_amd64.deb";
sudo dpkg -c cuda-repo-ubuntu1204_${CUDA_VERSION}_amd64.deb;
echo hi
sudo dpkg -i cuda-repo-ubuntu1204_${CUDA_VERSION}_amd64.deb;
find ./etc/apt/sources.list.d/cuda.list
sudo apt-get update -qq;
export CUDA_APT=${CUDA_VERSION%-*};
export CUDA_APT=${CUDA_APT/./-};
#sudo apt-get install -y cuda-drivers cuda-core-${CUDA_APT} cuda-cudart-dev-${CUDA_APT} cuda-cufft-dev-${CUDA_APT};
sudo apt-get install -y cuda;
sudo apt-get clean;
export CUDA_HOME=/usr/local/cuda-${CUDA_VERSION%%-*};
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH};
export PATH=${CUDA_HOME}/bin:${PATH};

