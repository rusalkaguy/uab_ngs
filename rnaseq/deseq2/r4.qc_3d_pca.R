# rgl interactive 3D
# or you could use the rgl package if you want to be able to rotate the plot.
# rgl example: https://support.bioconductor.org/p/33758/

# source("http://bioconductor.org/biocLite.R")  # get bioconductor
# biocLite("rgl")
library(rgl)

#Warning message:
#package .rgl. was built under R version 3.2.3 
rgl.open()
offset <- 50
par3d(windowRect=c(offset, offset, 640+offset, 640+offset))
rm(offset)
rgl.clear()
rgl.viewpoint(theta=45, phi=30, fov=60, zoom=1)
spheres3d(pc3$x[,1], pc3$x[,2], pc3$x[,3], radius=1.0, alpha=1, shininess=20)
aspect3d(1, 1, 1)
axes3d(col='black')
title3d("", "", "PC1", "PC2", "PC3", col='black')
#bg3d("white")
bg3d("gray")
rgl.clear(type='lights')
rgl.light(-45, 20, ambient='black', diffuse='#dddddd', specular='white')
rgl.light(60, 30, ambient='#dddddd', diffuse='#dddddd', specular='black')
rgl.texts(pc3$x[,1], pc3$x[,2], pc3$x[,3], colData(rld)$cond <- abbrev)
