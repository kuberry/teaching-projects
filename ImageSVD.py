import numpy
from scipy import linalg
import matplotlib
import matplotlib.pyplot as plt
import PIL
import Image

tolerance = 1e-3;

FILEIN='BMurray.jpg' #image can be in gif jpeg or png format 
FILEOUT='BMurray_out.jpg'
im=Image.open(FILEIN).convert('RGB')
pix=im.load()
w=im.size[0]
h=im.size[1]

M = numpy.zeros(shape=(w,h,3))
for k in range(3):
    for i in range(w):
        for j in range(h):
            M[i,j,k] = pix[i,j][k]

    U, s, Vh = linalg.svd(M[:,:,k], full_matrices=False)

    print U.shape, s.shape, Vh.shape, w, h
    s_total = s.sum()
    s /= s_total
    for key, val in enumerate(s):
        if val < tolerance:
            s[key] = 0;
    s *= s_total
    S = linalg.diagsvd(s, w, w)
    M[:,:,k] = numpy.dot(U, numpy.dot(S, Vh))

    
for i in range(w):
    for j in range(h):
        pix[i,j] =  (int(numpy.round(M[i,j,0])),int(numpy.round(M[i,j,1])),int(numpy.round(M[i,j,2])))

im.save(FILEOUT)






    

# Trigger emacs to run this script using the "compilye" command
# ;;; Local Variables: ***
# ;;; compile-command: "python ImageSVD.py" ***
# ;;; end: ***
