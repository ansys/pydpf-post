
@echo on

rem set ANSYS_PATH=C:\Program Files\ANSYS Inc\v211
rem set DPF_CONFIGURATION=release
rem set DPF_PATH=C:\Program Files\ANSYS Inc\v211\aisol\bin\winx64\Ans.Dpf.Grpc.exe
rem set ANS_PROTOCOL_ROOT=D:\AnsysDev\Protocols
set DPF_CORE_PATH=D:\AnsysDev\DPF-Core
set DPF_POST_PATH=D:\AnsysDev\DPF-Post

rem %ANS_PROTOCOL_ROOT%\packages\python\dpf;

set PYTHONPATH=%DPF_POST_PATH%

set root=C:\ProgramData\Anaconda3
call %root%\Scripts\activate.bat %root%

call jupyter lab

 

