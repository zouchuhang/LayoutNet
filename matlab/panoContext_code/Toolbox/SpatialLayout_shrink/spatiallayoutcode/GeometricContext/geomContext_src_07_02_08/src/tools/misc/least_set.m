function [ ptab, b, c, d, eps, ierror ] =  least_set ( ntab, xtab, ytab, ndeg )

%% LEAST_SET constructs the least squares polynomial approximation to data.
%
%  Discussion:
%
%    The least squares polynomial is not returned directly as a simple
%    polynomial.  Instead, it is represented in terms of a set of
%    orthogonal polynomials appopriate for the given data.  This makes
%    the computation more accurate, but means that the user can not
%    easily evaluate the computed polynomial.  Instead, the routine 
%    LEAST_EVAL should be used to evaluate the least squares polynomial
%    at any point.  (However, the value of the least squares polynomial
%    at each of the data points is returned as part of this computation.)
%
%
%    A discrete unweighted inner product is used, so that
%
%      ( F(X), G(X) ) = sum ( 1 <= I <= NTAB ) F(XTAB(I)) * G(XTAB(I)).
%
%    The least squares polynomial is determined using a set of
%    orthogonal polynomials PHI.  These polynomials can be defined
%    recursively by:
%
%      PHI(0)(X) = 1
%      PHI(1)(X) = X - B(1)
%      PHI(I)(X) = ( X - B(I) ) * PHI(I-1)(X) - D(I) * PHI(I-2)(X)
%
%    The array B(1:NDEG) contains the values
%
%      B(I) = ( X*PHI(I-1), PHI(I-1) ) / ( PHI(I-1), PHI(I-1) )
%
%    The array D(2:NDEG) contains the values
%
%      D(I) = ( PHI(I-1), PHI(I-1) ) / ( PHI(I-2), PHI(I-2) )
%
%    Using this basis, the least squares polynomial can be represented as
%
%      P(X)(I) = sum ( 0 <= I <= NDEG ) C(I) * PHI(I)(X)
%
%    The array C(0:NDEG) contains the values
%
%      C(I) = ( YTAB(I), PHI(I) ) / ( PHI(I), PHI(I) )
%
%  Modified:
%
%    16 May 2004
%
%  Reference:
%
%    Gisela Engeln-Muellges and Frank Uhlig,
%    Numerical Algorithms with C, pages 191-193.
%    Springer, 1996.
%
%  Parameters:
%
%    Input, integer NTAB, the number of data points.
%
%    Input, real XTAB(NTAB), the X data.  The values in XTAB
%    should be distinct, and in increasing order.
%
%    Input, real YTAB(NTAB), the Y data values corresponding
%    to the X data in XTAB.
%
%    Input, integer NDEG, the degree of the polynomial which the
%    program is to use.  NDEG must be at least 1, and less than or 
%    equal to NTAB-1.
%
%    Output, real PTAB(NTAB), the value of the least squares polynomial 
%    at the points XTAB(1:NTAB).
%
%    Output, real B(1:NDEG), C(1:NDEG+1), D(1:NDEG-1), arrays containing 
%    data about the polynomial.
%
%    Output, real EPS, the root-mean-square discrepancy of the
%    polynomial fit.
%
%    Output, integer IERROR, error flag.
%    zero, no error occurred;
%    nonzero, an error occurred, and the polynomial could not be computed.
%
  ierror = 0;
  C_OFFSET = 1;
  D_OFFSET = -1;
%
%  Check NDEG.
%
  if ( ndeg < 1 )
    ierror = 1;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'LEAST_SET - Fatal error!\n' );
    fprintf ( 1, '  NDEG < 1.\n' );
    error ( 'LEAST_SET - Fatal error!' );
  end

  if ( ntab <= ndeg )
    ierror = 1;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'LEAST_SET - Fatal error!\n' );
    fprintf ( 1, '  NTAB <= NDEG.\n' );
    error ( 'LEAST_SET - Fatal error!' );
  end
%
%  Check that the abscissas are strictly increasing.
%
  for ( i = 1 : ntab-1 )
    if ( xtab(i+1) <= xtab(i) )
      ierror = 1;
      fprintf ( 1, '\n' );
      fprintf ( 1, 'LEAST_SET - Fatal error!\n' );
      fprintf ( 1, '  XTAB must be strictly increasing, but\n' );
      fprintf ( 1, '  XTAB(%d) = %f\n', i,   xtab(i)   );
      fprintf ( 1, '  XTAB(%d) = %f\n', i+1, xtab(i+1) );
      error ( 'LEAST_SET - Fatal error!' );
    end
  end

  i0l1 = 0;
  i1l1 = ntab;
%
%  The polynomial is of degree at least 0.
%
  y_sum = sum ( ytab(1:ntab) );
  rn0 = ntab;
  c(0+C_OFFSET) = y_sum / ntab;

  ptab(1:ntab) = y_sum / ntab;

  if ( ndeg == 0 )
    eps = sum ( ( ptab(1:ntab) - ytab(1:ntab) ).^2 );
    eps = sqrt ( eps / ntab );
    b = [];
    d = [];
    return;
  end
%
%  The polynomial is of degree at least 1.
%
  b(1) = sum ( xtab(1:ntab) ) / ntab;

  s = 0.0E+00;
  sum2 = 0.0E+00;
  for ( i = 1 : ntab )
    ztab(i1l1+i) = xtab(i) - b(1);
    s = s + ztab(i1l1+i)^2;
    sum2 = sum2 + ztab(i1l1+i) * ( ytab(i) - ptab(i) );
  end

  rn1 = s;
  c(1+C_OFFSET) = sum2 / s;

  for ( i = 1 : ntab )
    ptab(i) = ptab(i) + c(1+C_OFFSET) * ztab(i1l1+i);
  end

  if ( ndeg == 1 )
    eps = sum ( ( ptab(1:ntab) - ytab(1:ntab) ).^2 );
    eps = sqrt ( eps / ntab );
    d = [];
    return;
  end

  ztable(1:ntab) = 1.0E+00;

  mdeg = 2;
  k = 2;

  while ( 1 )

    d(k+D_OFFSET) = rn1 / rn0;

    sum2 = 0.0E+00;
    for ( i = 1 : ntab )
      sum2 = sum2 + xtab(i) * ztab(i1l1+i) * ztab(i1l1+i);
    end

    b(k) = sum2 / rn1;

    s = 0.0E+00;
    sum2 = 0.0E+00;
    for ( i = 1 : ntab )
      ztab(i0l1+i) = ( xtab(i) - b(k) ) * ztab(i1l1+i) ...
        - d(k+D_OFFSET) * ztab(i0l1+i);
      s = s + ztab(i0l1+i) * ztab(i0l1+i);
      sum2 = sum2 + ztab(i0l1+i) * ( ytab(i) - ptab(i) );
    end

    rn0 = rn1;
    rn1 = s;
 
    c(k+C_OFFSET) = sum2 / rn1;

    it = i0l1;
    i0l1 = i1l1;
    i1l1 = it;

    for ( i = 1 : ntab )
      ptab(i) = ptab(i) + c(k+C_OFFSET) * ztab(i1l1+i);
    end
    
    if ( ndeg <= mdeg )
      break;
    end

    mdeg = mdeg + 1;
    k = k + 1;

  end
%
%  Compute the RMS error.
%
  eps = sum ( ( ptab(1:ntab) - ytab(1:ntab) ).^2 );
  eps = sqrt ( eps / ntab );
