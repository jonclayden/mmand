# Mathematical Morphology in Any Number of Directions

The `mmand` R package provides tools for performing mathematical morphology operations, such as erosion and dilation, on arrays of arbitrary dimensionality. It can also resample arrays, obtaining values between pixel centres or scaling the image up or down wholesale.

A test image of a jet engine fan is available within the package, and will be used for demonstration below. It can be read in and displayed using the code

```R
library(mmand)
library(png)

fan <- readPNG(system.file("images", "fan.png", package="mmand"))
display(fan)
```

![Fan test image](http://www.clayden.org/files/mmand/fan.png)

## Mathematical morphology

[Mathematical morphology](http://en.wikipedia.org/wiki/Mathematical_morphology) is an image processing technique that can be used to emphasise or remove certain types of features from binary or greyscale images. It is classically performed on two-dimensional images, but can be useful in three or more dimensions, as, for example, in medical image analysis. The `mmand` package can work on R arrays of any dimensionality, including one-dimensional vectors.

The basic operations in mathematical morphology are *erosion* and *dilation*. A simple one-dimensional example serves to illustrate their effects:

```R
x <- c(0,0,1,0,0,0,1,1,1,0,0)
k <- c(1,1,1)

erode(x,k)
# [1] 0 0 0 0 0 0 0 1 0 0 0

dilate(x,k)
# [1] 0 1 1 1 0 1 1 1 1 1 0
```

The `erode()` function "thins out" areas in the input vector, `x`, which were "on" (i.e. set to 1), to the point where the first of these areas "disappears" entirely. Conversely, the `dilate()` function expands these regions into neighbouring pixels.

The vector `k` here is called the *kernel* or *structuring element*. It effectively controls the region of influence of the operation when it is applied to each value.

Derived from these basic operations are the *opening* and *closing* functions. These apply both basic operations, using the same kernel, but in different orders: the opening is an erosion followed by a dilation, whereas a closing is a dilation followed by an opening.

```R
opening(x,k)
# [1] 0 0 0 0 0 0 1 1 1 0 0

closing(x,k)
# [1] 0 0 1 0 0 0 1 1 1 0 0
```

Notice that, in this case, the closing gets us back to where we started, whereas the opening does not. This is because the initial erosion operation removes the first "on" block entirely, so it cannot be recovered by the subsequent dilation.

## Greyscale morphology

Mathematical morphology is not limited to binary data. When generalised to greyscale images, erosion replaces each nonzero pixel with the minimum value within the kernel when it is centred at that pixel, and dilation uses the maximum. For example,

```R
x <- c(0,0,0.5,0,0,0,0.2,0.5,0.3,0,0)
erode(x,k)
# [1] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.0 0.0 0.0
```

Notice that the remaining nonzero value is now reduced from 0.5 to 0.2, the minimum value across the original pixel and its neighbours on either side. With a wider kernel, its final value would have dropped to zero.
