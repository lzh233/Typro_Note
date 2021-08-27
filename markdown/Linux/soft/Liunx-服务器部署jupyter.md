# 服务器部署jupyter
0. Add conda path
```shell
 export PATH="/SGRNJ03/randd/public/soft/conda/miniconda/bin:$PATH" 
 source ~/.bashrc`
```

```shell
#安装jupyter
# Assuming your conda-env is named cenv
conda create -n cenv
conda activate cenv
conda install -c anaconda jupyter
```

1. Install python kernel
``` shell
conda install ipykernel
ipython kernel install --user --name=<any_name_for_kernel>
```
3. Install R kernel
```shell
conda install -c r r-irkernel
IRkernel::installspec(name="<any_name_for_kernel>",="<any_name_for_kernel>")
```
4.Run

``` shell
mkdir ~/.jupyter/
cp /SGRNJ/Database/script/pipe/develop/shared/config/jupyter_notebook_config.py ~/.jupyter/
mkdir /SGRNJ03/randd/user/{user}/jupyter && cd /SGRNJ03/randd/user/{user}/jupyter
nohup jupyter notebook &
```

5.updateR

```shell
conda install -c r r-base=4.0.3
#or
conda install R=4.0
```

