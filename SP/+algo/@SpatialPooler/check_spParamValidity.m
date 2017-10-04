function check_spParamValidity(spParams,ncolumns,ninput)

assert(ncolumns>=sum(ninput), 'SP to input mapping doesn''t work if input size is larger than SP size');
