clc;


pth=["/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/D28";           %1
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/D50";            %2
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/D200";           %3
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/D100";           %4
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/U50";            %5
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/U100";           %6
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/U200";           %7
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/Q50-100";        %8
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/Q100-100";       %9
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/Q100-150";       %10
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/Q100-200";       %11
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Thesis/Q200-100";       %12
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/A-D40.8";    %13
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/A-D46"       %14
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/A-D54";      %15
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/A-Q100-150"; %16
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/A-Q100-200"; %17
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/B-D52.2";    %18
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/B-Q100-150"; %19
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/C-D47";      %20
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/C-Q100-150"; %21    
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/A-D28";      %22
    "/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Simulations_Article_2/A-Q100-100"];%23    


tgPath="/media/geomas/86e0a188-ce71-4076-82d4-812cf66dd378/Simulations/Res";
mode=[ "LOOPVR","LOOPS","FORCEXT","VOID"];
%"ANI","EDG","FC","FCCL","FCCLTF","FORCE","FORCEXT","LOOPS","LOOPA",
%    "LOOPCT","LOOPD","LOOPVR","VOID","STRAIN","STNW2AC"


for i=23:23
    exe_NoApp(pth(i))
    resGrouper(tgPath,pth(i),mode)
end

%quit

%Check if run on terminal, if Yes quit matlab
%if ~mclIsNoDisplaySet()
%    return
%end

