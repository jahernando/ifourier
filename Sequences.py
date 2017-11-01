#
# Module with useful functions for secuences.
#
#
# Author: Gonzalo Martinez
#
# Last update: 11 / 04 / 2014
#

from math import ceil
from copy import deepcopy

def Zeros( rows, cols = 0 ):
    ''' This function makes a vector (matrix) of size rows ( rows x cols ) filled with zeros.'''
    
    return [ 0. for i in range(rows) ] if not cols else [ Zeros(cols) for i in range(rows) ]

def Ones( rows, cols=0 ):
    ''' This function makes a vector (matrix) of size rows ( rows x cols ) filled with ones.'''
    
    return [ 1. for i in range(rows) ] if not cols else [ Ones(cols) for i in range(rows) ]

def Identity( N ):
    ''' This functions returns an identity matrix of size N x N.'''

    return [ [1. if i==j else 0. for j in xrange(N)] for i in xrange(N) ]

def Scale( M, F ):
    ''' This function takes a vector or a matrix of numbers and re-scale it by a factor F.'''
    
    return map( lambda x: Scale( x, F ), M ) if isinstance( M[0], (list,tuple) ) else map( lambda y: y * F, M )

def Multiply( A, B ):
    ''' Multiplies element by element two sequences.'''
    
    if not isinstance( A, list ):
        return A * B

    if not len(A) == len(B):
        raise ValueError('Arguments must have the same size')

    C = deepcopy(A)
    for i in range(len(B)):
        C[i] = Multiply( C[i], B[i] )

    return C

def Apply( A, F ):
    ''' Applies the function F to all the elements of A.'''

    if not isinstance( A, list ):
        return F( A )
    
    B = deepcopy( A )
    
    for i in range(len(B)):
        B[i] = Apply( B[i], F )

    return B

def applyto(A,F):
    ns = A.shape
    B = deepcopy(A)
    if (len(ns)==1):
        for i in range(ns[0]):
            B[i]=F(A[i])
    if (len(ns)==2):
        for i in range(ns[0]):
            for j in range(ns[1]):
                B[i][j]=F(A[i][j])
    if (len(ns)==3):
        for i in range(ns[0]):
            for j in range(ns[1]):
                for k in range(ns[2]):
                    B[i][j][k]=F(A[i][j][k])
    return B

def PrintMatrix( M ):
    ''' This function prints a matrix in a more visual form.'''
    
    for line in M:
        print line
    
    return None

def Transpose(M):
    MT = deepcopy(M)
    n,m = len(M),len(M[0])
    for i in range(n):
        for j in range(m):
            MT[j][i]=M[i][j]
    return MT

def Swap( L, i, j ):
    ''' Swap two values in a list.'''
    
    L[i],L[j] = L[j],L[i]

def Sort( l, l2 = None ):
    ''' Sorts a list in ascending order.'''
    
    if isinstance( l2, list ):
        return Sort2( l, l2 )
    
    elif isinstance( l2, (int,float) ):
        return SortNumber( l, l2 )
    
    
    L = list( l )
    for i in range( len(L) - 1 ):
        A = L[ i + 1 ]
        j = i
        
        while j+1 and L[j] > A:
            L[ j+1 ] = L[ j ]
            j -= 1

        L[ j+1 ] = A

    return L

def Sort2( l1, l2 ):
    ''' Sorts the first list in ascending order and the second one with the order of the former one.'''
    
    L1 = list( l1 )
    L2 = list( l2 )
    
    for i in range( len(L1) - 1 ):
        A = L1[ i + 1 ]
        B = L2[ i + 1 ]
        j = i
        
        while j+1 and L[j] > A:
            L1[ j+1 ] = L1[ j ]
            L2[ j+1 ] = L2[ j ]
            j -= 1
        
        L1[ j+1 ] = A
        L2[ j+1 ] = B
    
    return L1, L2
    
def SortNumber( L, number = 0 ):
    ''' This function orders a list by its distance to a given number.'''
    
    return list( zip( *sorted( zip( map( lambda x: abs( x - number ), L ), L ) ) )[1])

def Reversed( L ):
    ''' This function reverses the order of a list without destroying the original.'''
    
    return L[::-1]

def Izip( Z ):
    ''' This function inverts the zipping in a list of pairs without destroying the original.'''
    
    return zip( Reversed( zip( *Z ) ) )

def MaxIndex( x ):
    ''' This function returns the indices of the greatest value in a list.'''
    
    if isinstance( x[0], (int,float) ):
        maxval = max(x)
        return x.index( maxval), maxval
    
    indices, values = zip( *map( MaxIndex, x ) )
    i0, maxval = MaxIndex( values )
    i1 = indices[ i0 ]
    return i0, i1, maxval

def MinIndex( x ):
    ''' This function returns the indeices of the smallest value in a list.'''
    
    if isinstance( x[0], (int,float) ):
        minval = min(x)
        return x.index( minval), minval
    
    indices, values = zip( *map( MinIndex, x ) )
    i0, minval = MinIndex( values )
    i1 = indices[ i0 ]
    return i0, i1, minval

def Cumulative( L ):
    ''' This function return the cumulative of a given list.'''
    
    R = [0]
    for l in L:
        R.append( R[-1] + l )

    return R[1:]

def Frecs( l ):
    ''' This function takes a list and returns a dictionary with the frecuency of each element.'''
    
    frecs={}
    for element in l:
        if element in frecs:
            continue
        frecs[element] = l.count( element )
    
    return frecs


def Binning( nbins = 100, lower = 0., upper = 1.):
    ''' Generalization of range: returns a inverval sliced in several partitions (bins). The return is a list which values are the beginning of each bin. The arguments are the number of bins and the minimum and the maximum of the interval.'''
    
    nbins = int( ceil( nbins ) )
    size = float( upper - lower )/ nbins
    
    return [ lower + i * size for i in range( int(nbins) ) ]


def mapA2(f,A):
    return f(map(f,A))

def mapA3(f,A):
    return f(map(lambda ai: mapA2(f,ai),A))
