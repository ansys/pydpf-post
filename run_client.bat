
@echo on

set DPF_POST_PATH=D:\AnsysDev\DPF-Post
set PYTHONPATH=%DPF_POST_PATH%

set root=C:\ProgramData\Anaconda3
call %root%\Scripts\activate.bat %root%

call jupyter lab

 

