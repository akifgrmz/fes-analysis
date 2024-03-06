function [OutSig,ClippedSig]=amp_modul(InSig,scale,offset,asf)
% Adaptive amplitude modulation for low pass filtering MAV signals

def_asf=[1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,5,5,6,7,8,9,9,10,11,11];
% def_asf=[1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4,4,5,5,6,7,8,9,9,10,11,11];
% def_asf=[1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,5,5,6,7,8,9,9,10,11,11];
% def_asf=[1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,8,10];
% def_asf=[1,2,3,4,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1];
% def_asf=[1 3 7 15 31 35 44 38 29 21 15 8 5 3 2 1 1 1 1];
% def_asf=[1 3 7 15 25 33 38 33 25 18 13 8 5 3 2 1 1];

def_offset=1;
def_scale=0;
dir=0;
asf_index=1;

if nargin == 0
    error('Not enough inputs! ')
elseif nargin == 1 
    asf=def_asf;
    scale=def_scale;
    offset=def_offset;
elseif nargin == 2 
    asf=def_asf;
    scale=def_scale;
elseif nargin == 3 
    asf=def_asf;
elseif nargin >4 
    error('Too many inputs! ')
end
 
max_index=length(asf);
ClippedSig=clip((InSig-offset)*scale,0,1);
OutSig=ClippedSig*0;
OutSig(1)=ClippedSig(1);

for iSamp=2:length(ClippedSig)
    delta=ClippedSig(iSamp)-OutSig(iSamp-1);
    new_dir=sign(delta);
    mag=abs(delta);

    if new_dir==0
        asf_index=1;
    elseif new_dir==dir
        asf_index=asf_index+1;

        if asf_index>max_index
            asf_index=max_index;
        end
    else
        asf_index=1;
    end
    
    dir=new_dir;
    step=min(mag,asf(asf_index))*dir;
    OutSig(iSamp)=OutSig(iSamp-1)+step;
    
    OutSig=clip(OutSig,0,0.15);
end

function[clipped]=clip(input,low,high)
    clipped=input*0;
    for i=1:length(input)
        clipped(i)=min(max(input(i),low),high);
    end
end

end