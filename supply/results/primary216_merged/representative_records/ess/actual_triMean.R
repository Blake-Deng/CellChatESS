function (x, na.rm = TRUE) 
{
    mean(collapse::fquantile(x, probs = c(0.25, 0.5, 0.5, 0.75), 
        na.rm = na.rm))
}
