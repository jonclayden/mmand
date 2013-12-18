boxKernel <- function ()
{
    return (structure(list(name="box"), class="kernel"))
}

triangleKernel <- function ()
{
    return (structure(list(name="triangle"), class="kernel"))
}

mitchellNetravaliKernel <- mnKernel <- function (B = 1/3, C = 1/3)
{
    return (structure(list(name="mitchell-netravali", B=B, C=C), class="kernel"))
}

sampleKernel <- function (kernel, values)
{
    if (!is.list(kernel) || !("kernel" %in% class(kernel)))
        output(OL$Error, "Specified kernel is invalid")
    
    return (.Call("sample_kernel", kernel, as.numeric(values), PACKAGE="irk"))
}
