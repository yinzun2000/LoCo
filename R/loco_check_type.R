loco_check_type <-  function(x){
    type.o  <-  typeof(x)
    if (type.o == 'integer' |
        type.o == 'double')
    {
        type.o  <-  'n'
    } else if (type.o == 'character')
    {
        type.o  <-  'c'
    } else if (type.o == 'logical')
    {
        type.o  <-  'l'
    } else
        type.o  <-  'o'

    dim.o   <-  dim(x)
    if (is.null(dim.o))
        dim.o   <-  length(x)

    v.o <-  list(t=type.o,d=dim.o)
    return(v.o)
}
