import numpy as np

# k x0 xf y0 yf
# 0 0 300 0 300
# 1 0 300 210 510
# 2 0 300 420 720
# 3 0 300 630 930
# 4 0 300 840 1140
# 5 210 510 0 300
# 6 210 510 210 510
# 7 210 510 420 720
# 8 210 510 630 930
# 9 210 510 840 1140
# 10 420 720 0 300
# 11 420 720 210 510
# 12 420 720 420 720
# 13 420 720 630 930
# 14 420 720 840 1140
# 15 630 930 0 300
# 16 630 930 210 510
# 17 630 930 420 720
# 18 630 930 630 930
# 19 630 930 840 1140
# 20 840 1140 0 300
# 21 840 1140 210 510
# 22 840 1140 420 720
# 23 840 1140 630 930
# 24 840 1140 840 1140
def scanning(data,Lsiz=128,cut=29):
    """
    This function returns an array with the subpatches of the focused
    and defocused images
    """
    #Lsiz=int(N) #Subframe size
    lsiz=int(Lsiz-2*cut) #Size of the subframe that is not overlapped to others
    if np.floor(data.shape[0]/lsiz)>2:
        i_max=int(np.floor(data.shape[0]/lsiz)-1)
    else:#To perform a 2 x 2 subfielding
        i_max=int(np.floor(data.shape[0]/lsiz))
    i_vec=np.arange(0,i_max)

    #Overlapping of subframes
    kk=-1
    if data.ndim<3:
        Nima=1
        data_subpatch=np.zeros((int(i_max**2),Lsiz,Lsiz))
    else:    
        Nima=data.shape[2]
        data_subpatch=np.zeros((int(i_max**2),Lsiz,Lsiz,Nima))
    for i in i_vec:
        xc=int(Lsiz/2+lsiz*i)
        x0=int(xc-Lsiz/2)
        xf=int(xc+Lsiz/2)
        for j in i_vec:
            kk+=1
            yc=int(Lsiz/2+lsiz*j)
            y0=int(yc-Lsiz/2)
            yf=int(yc+Lsiz/2)
            print(kk,x0,xf,y0,yf)
            if x0<0 or y0<0:
                print('Error in scanning.py: x0 or y0 cannot be negative')
                quit()
            elif xf>data.shape[0] or yf>data.shape[1]:
                print('Error in scanning.py: xf or yf larger than of size')
                print('xf:',xf)
                print('yf:',yf)
            if Nima==1:
                data_subpatch[kk,:,:]=data[x0:xf,y0:yf]    
            else:
                for n in range(Nima):
                    data_subpatch[kk,:,:,n]=data[x0:xf,y0:yf,n]
    return data_subpatch