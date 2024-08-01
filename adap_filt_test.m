function [out,ClippedSig,dir,asf_index]=adap_filt_test(InSig,scale,offset,asf,dir,asf_index)
% Adaptive amplitude modulation for low pass filtering MAV signals

def_asf=[1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,5,5,6,7,8,9,9,10,11,11];


def_offset=1;
def_scale=0;

max_index=length(asf);
ClippedSig=clip((InSig-offset)*scale,0,1);
OutSig=ClippedSig*0;
OutSig(1)=ClippedSig(1);

    delta=ClippedSig(2)-OutSig(1);
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
    OutSig(2)=OutSig(1)+step;
    
    OutSig=clip(OutSig,0,.15);
out=OutSig(2);

function[clipped]=clip(input,low,high)
    clipped=input*0;
    for i=1:length(input)
        clipped(i)=min(max(input(i),low),high);
    end
end

end