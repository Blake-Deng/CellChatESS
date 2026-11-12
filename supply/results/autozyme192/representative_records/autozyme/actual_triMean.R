function (...) 
{
    if (is_disabled()) {
        return(.zyme_strip_and_forward(original, ...))
    }
    fast_fn(...)
}
