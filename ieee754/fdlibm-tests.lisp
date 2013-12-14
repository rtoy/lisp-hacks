(in-package "FDLIBM")

(rt:deftest sin.1
    (%sin 0d0)
  0d0)

(rt:deftest sin.2
    (%sin -0d0)
  -0d0)

(rt:deftest sin.3
    ;; Tests the case for |x| < 2^-27, but not 0.
    (%sin (scale-float 1d0 -28))
  #.(scale-float 1d0 -28))

(rt:deftest sin.4
    ;; Just a random test, without argument reduction
    (%sin .5d0)
  0.479425538604203d0)

(rt:deftest sin.5
    ;; Test for arg near pi/2
    (%sin (/ pi 2))
  1d0)

(rt:deftest sin.red.0
    ;; Test for argument reduction with n mod 4 = 0
    (%sin (* 7/4 pi))
  -7.07106781186547675943154203316156531867416581156d-1)

(rt:deftest sin.red.1
    ;; Test for argument reduction with n mod 4 = 1
    (%sin (* 9/4 pi))
  7.07106781186547329560731709118834541043171055432d-1)

(rt:deftest sin.red.2
    ;; Test for argument reduction with n mod 4 = 2
    (%sin (* 11/4 pi))
  7.07106781186548390575743300374993861263439430213d-1)

(rt:deftest sin.red.3
    ;; Test for argument reduction with n mod 4 = 3
    (%sin (* 13/4 pi))
  -7.07106781186547871002109559079472349116005337743d-1)

(rt:deftest sin.misc.1
    ;; Test for argument reduction
    (%sin (scale-float 1d0 120))
  0.377820109360752d0)

(rt:deftest cos.1
    (%cos 0d0)
  1d0)

(rt:deftest cos.2
    (%cos -0d0)
  1d0)

(rt:deftest cos.3
    ;; Test for |x| < 2^-27
    (%cos (scale-float 1d0 -28))
  1d0)

(rt:deftest cos.6
    ;; Test for branch |x| < .3
    (%cos 0.25d0)
  0.9689124217106447d0)

(rt:deftest cos.5
    ;; Test for branch |x| > .3 and \x| < .78125
    (%cos 0.5d0)
  8.7758256189037271611628158260382965199164519711d-1)

(rt:deftest cos.6
    ;; Test for branch |x| > .3 and |x| > .78125
    (%cos 0.785d0)
  0.7073882691671998d0)

(rt:deftest cos.7
    ;; Random test near pi/2
    (%cos (/ pi 2))
  6.123233995736766d-17)

(rt:deftest cos.8
    ;; Test for argument reduction
    (%cos (scale-float 1d0 120))
  -0.9258790228548379d0)


