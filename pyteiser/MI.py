{
    const double *data;
    const int *nrows, *ncols, *nbins;
    int *res;
    double *spl, *col;
    SEXP Rres, Rspl, Rcol;
    PROTECT(Rdata=AS_NUMERIC(Rdata));
    PROTECT(Rnrows=AS_INTEGER(Rnrows));
    PROTECT(Rncols=AS_INTEGER(Rncols));
    PROTECT(Rnbins=AS_INTEGER(Rnbins));
    data = NUMERIC_POINTER(Rdata);
    nrows = INTEGER_POINTER(Rnrows);
    ncols = INTEGER_POINTER(Rncols);
    nbins = INTEGER_POINTER(Rnbins);
    PROTECT(Rres=NEW_INTEGER((*nrows)*(*ncols)));
    PROTECT(Rspl=NEW_NUMERIC(*nbins));
    PROTECT(Rcol=NEW_NUMERIC(*nrows));
    spl = NUMERIC_POINTER(Rspl);
    col = NUMERIC_POINTER(Rcol);
    res = INTEGER_POINTER(Rres);
    for( int v=0; v<(*ncols); ++v )
    {
        int N=(*nrows);
        for( int s=0; s<N; ++s )
            col[s] = data[v*N+s];
            std::sort(col,col+N);
		//for( int j=0; ISNA(col[j]); ++j ) N--;
        for( int j=N-1; j > 0 && ISNA(col[j]); --j ) N--;
        int freq = N/(*nbins), mod = N%(*nbins);
        int splitpoint=freq-1;
        for( int i=0; i<(*nbins)-1; ++i ) {
              if( mod>0 ) {spl[i] = col[splitpoint+1]; mod--;}
              else spl[i]=col[splitpoint];
              splitpoint += freq;
         }
         spl[(*nbins)-1] = col[N-1]+epsilon;
         for( int s=0; s<(*nrows); ++s )
            if( !ISNA(data[s+v*(*nrows)]) )
            {
                int bin = -1;
                for( int k=0; bin==-1 && k<(*nbins); ++k )
                   if( data[s+v*(*nrows)] <= spl[k] ) bin=k;
                res[s+v*(*nrows)] = bin+1;
            }
            else res[s+v*(*nrows)]= NA_INTEGER;
     }
     UNPROTECT(7);
     return Rres;
}

def discretize_eq_freq(data, nbins):

    nrows = data.shape[0]
    ncols = data.shape[1]

    col = []

    res = np.zeros(data.shape)
    EPSILON = 0.01

    for v in range(ncols):
        N = nrows
        for s in range(nrows):
            col[s] = data[v*N + s]
