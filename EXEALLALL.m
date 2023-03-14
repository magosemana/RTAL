%clc
%{
pth=["~/Documents/Joao/3D/3P2022-03-09/D28";    %1
    "~/Documents/Joao/3D/3P2022-03-09/D50";     %2
    "~/Documents/Joao/3D/3P2022-03-09/D100";    %3
    "~/Documents/Joao/3D/3P2022-03-09/D200";    %4
    "~/Documents/Joao/3D/3P2022-03-09/U50";     %5
    "~/Documents/Joao/3D/3P2022-03-09/U100";    %6
    "~/Documents/Joao/3D/3P2022-03-09/U200";    %7
    "~/Documents/Joao/3D/3P2022-08-10/Q50";     %8
    "~/Documents/Joao/3D/3P2022-08-10/Q100";    %9
    "~/Documents/Joao/3D/3P2022-08-10/Q200";    %10
    "~/Documents/Joao/3D/3P2022-08-10/Q100-200";%11
    "~/Documents/Joao/3D/3P2022-08-10/Q100-150";%12
    "~/Documents/Joao/3D/3P2022-12-22/D100";    %13
    "~/Documents/Joao/3D/3P2023-02-14/Q50-100"; %14
    "~/Documents/Joao/3D/3P2023-02-14/Q100-100"; %15
    "~/Documents/Joao/3D/3P2023-02-14/Q100-150"; %16
    "~/Documents/Joao/3D/3P2023-02-14/Q100-200"; %17
    "~/Documents/Joao/3D/3P2023-02-14/Q200-100"];%18
%}

pth=["/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/D28";      %1
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/D50";       %2
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/D100";      %3
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/D200";      %4
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/U50";       %5
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/U100";      %6
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/U200";      %7
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Q100-100"]; %8 

tgPath="~/Documents/Joao/3D/Res";
mode=["STRAIN","FC","FCCL","LOOPS","LOOPA","LOOPD","LOOPVR","VOID"];
%"ANI","EDG","FC","FCCL","FCCLTF","FORCE","FORCEXT","LOOPS","LOOPA",
%    "LOOPCT","LOOPD","LOOPVR","VOID","STRAIN"


for i=1:3%numel(pth)
    exe_NoApp(pth(i))
    %resGrouper(tgPath,pth(i),mode)
end

%quit

%Check if run on terminal, if Yes quit matlab
%if ~mclIsNoDisplaySet()
%    return
%end

