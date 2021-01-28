% generates the mex executable
GTRACK_PATH = 'C:\Radar\top\mmwave_sdk\ti\alg\gtrack';
GTRACK_INCLUDE_PATH = GTRACK_PATH;
GTRACK_LIB_PATH = [GTRACK_PATH,'\lib'];

%cmd = 'mex -D__GNUC__ -Dmex -DWIN32 -D_WIN32 ';
cmd = 'mex -g ';
cmd = [cmd ['-I"',GTRACK_INCLUDE_PATH,'" ']];
cmd = [cmd ['-L"',GTRACK_LIB_PATH,'" ']];
cmd = [cmd '-lgtrackMex '];
cmd = [cmd 'gtrack_create_mex.c '];
cmd = [cmd 'gtrack_osal.c '];
disp(cmd);
eval(cmd)

cmd = 'mex -g ';
cmd = [cmd ['-I"',GTRACK_INCLUDE_PATH,'" ']];
cmd = [cmd ['-L"',GTRACK_LIB_PATH,'" ']];
cmd = [cmd '-lgtrackMex '];
cmd = [cmd 'gtrack_step_mex.c '];
cmd = [cmd 'gtrack_osal.c '];
disp(cmd);
eval(cmd)


% mex '-LC:\Radar\gtrack\msvc\x64\Debug' '-IC:\Radar\gtrack\api' -lgtrackLib gtrack_moduleCreate_mex.c osal_memoryAlloc.c
