(in-package "FDLIBM")

(declaim (ext:start-block kernel-sin kernel-cos kernel-tan %sin %cos %tan))

;; kernel sin function on [-pi/4, pi/4], pi/4 ~ 0.7854
;; Input x is assumed to be bounded by ~pi/4 in magnitude.
;; Input y is the tail of x.
;; Input iy indicates whether y is 0. (if iy=0, y assume to be 0). 
;;
;; Algorithm
;;	1. Since sin(-x) = -sin(x), we need only to consider positive x. 
;;	2. if x < 2^-27 (hx<0x3e400000 0), return x with inexact if x!=0.
;;	3. sin(x) is approximated by a polynomial of degree 13 on
;;	   [0,pi/4]
;;		  	         3            13
;;	   	sin(x) ~ x + S1*x + ... + S6*x
;;	   where
;;	
;; 	|sin(x)         2     4     6     8     10     12  |     -58
;; 	|----- - (1+S1*x +S2*x +S3*x +S4*x +S5*x  +S6*x   )| <= 2
;; 	|  x 					           | 
;; 
;;	4. sin(x+y) = sin(x) + sin'(x')*y
;;		    ~ sin(x) + (1-x*x/2)*y
;;	   For better accuracy, let 
;;		     3      2      2      2      2
;;		r = x *(S2+x *(S3+x *(S4+x *(S5+x *S6))))
;;	   then                   3    2
;;		sin(x) = x + (S1*x + (x *(r-y/2)+y))

(declaim (ftype (function (double-float double-float fixnum)
			  double-float)
		kernel-sin))

(defun kernel-sin (x y iy)
  (declare (type (double-float -1d0 1d0) x y)
	   (fixnum iy)
	   (optimize (speed 3) (safety 0)))
  (let ((ix (ldb (byte 31 0) (kernel:double-float-high-bits x))))
    (when (< ix #x3e400000)
      (if (zerop (truncate x))
	  (return-from kernel-sin x)
	  (return-from kernel-sin x)))
    (let* ((s1 -1.66666666666666324348d-01)
	   (s2  8.33333333332248946124d-03)
	   (s3 -1.98412698298579493134d-04)
	   (s4  2.75573137070700676789d-06)
	   (s5 -2.50507602534068634195d-08)
	   (s6  1.58969099521155010221d-10)
	   (z (* x x))
	   (v (* z x))
	   (r (+ s2
		 (* z
		    (+ s3
		       (* z
			  (+ s4
			     (* z
				(+ s5
				   (* z s6))))))))))
      (if (zerop iy)
	  (+ x (* v (+ s1 (* z r))))
	  (- x (- (- (* z (- (* .5 y)
			     (* v r)))
		     y)
		  (* v s1)))))))

;; kernel cos function on [-pi/4, pi/4], pi/4 ~ 0.785398164
;; Input x is assumed to be bounded by ~pi/4 in magnitude.
;; Input y is the tail of x. 
;;
;; Algorithm
;;	1. Since cos(-x) = cos(x), we need only to consider positive x.
;;	2. if x < 2^-27 (hx<0x3e400000 0), return 1 with inexact if x!=0.
;;	3. cos(x) is approximated by a polynomial of degree 14 on
;;	   [0,pi/4]
;;		  	                 4            14
;;	   	cos(x) ~ 1 - x*x/2 + C1*x + ... + C6*x
;;	   where the remez error is
;;	
;; 	|              2     4     6     8     10    12     14 |     -58
;; 	|cos(x)-(1-.5*x +C1*x +C2*x +C3*x +C4*x +C5*x  +C6*x  )| <= 2
;; 	|    					               | 
;; 
;; 	               4     6     8     10    12     14 
;;	4. let r = C1*x +C2*x +C3*x +C4*x +C5*x  +C6*x  , then
;;	       cos(x) = 1 - x*x/2 + r
;;	   since cos(x+y) ~ cos(x) - sin(x)*y 
;;			  ~ cos(x) - x*y,
;;	   a correction term is necessary in cos(x) and hence
;;		cos(x+y) = 1 - (x*x/2 - (r - x*y))
;;	   For better accuracy when x > 0.3, let qx = |x|/4 with
;;	   the last 32 bits mask off, and if x > 0.78125, let qx = 0.28125.
;;	   Then
;;		cos(x+y) = (1-qx) - ((x*x/2-qx) - (r-x*y)).
;;	   Note that 1-qx and (x*x/2-qx) is EXACT here, and the
;;	   magnitude of the latter is at least a quarter of x*x/2,
;;	   thus, reducing the rounding error in the subtraction.
(declaim (ftype (function (double-float double-float)
			  double-float)
		kernel-cos))

(defun kernel-cos (x y)
  (declare (type (double-float -1d0 1d0) x y)
	   (optimize (speed 3) (safety 0)))
  ;; cos(-x) = cos(x), so we just compute cos(|x|).
  (let ((ix (ldb (byte 31 0) (kernel:double-float-high-bits x))))
    ;; cos(x) = 1 when |x| < 2^-27
    (when (< ix #x3e400000)
      ;; Signal inexact if x /= 0
      (if (zerop (truncate x))
	  (return-from kernel-cos 1d0)
	  (return-from kernel-cos 1d0)))
    (let* ((c1  4.16666666666666019037d-02)
	   (c2 -1.38888888888741095749d-03)
	   (c3  2.48015872894767294178d-05)
	   (c4 -2.75573143513906633035d-07)
	   (c5  2.08757232129817482790d-09)
	   (c6 -1.13596475577881948265d-11)
	   (z (* x x))
	   (r (* z
		 (+ c1
		    (* z
		       (+ c2
			  (* z
			     (+ c3
				(* z
				   (+ c4
				      (* z
					 (+ c5
					    (* z c6)))))))))))))
      (cond ((< ix #x3fd33333)
	     ;; \x| < 0.3
	     (- 1 (- (* .5 z)
		     (- (* z r)
			(* x y)))))
	    (t
	     (let* ((qx (if (> ix #x3fe90000)
			    0.28125d0
			    ;; x/4, exactly, and also dropping the
			    ;; least significant 32 bits of the
			    ;; fraction. (Why?)
			    (kernel:make-double-float (- ix #x00200000)
						      0)))
		    (hz (- (* 0.5 z) qx))
		    (a (- 1 qx)))
	       (- a (- hz (- (* z r)
			     (* x y))))))))))

(declaim (type (simple-array double-float (*)) ttt))
(defconstant ttt
  (make-array 13 :element-type 'double-float
	      :initial-contents
	      '(3.33333333333334091986d-01
		1.33333333333201242699d-01
		5.39682539762260521377d-02
		2.18694882948595424599d-02
		8.86323982359930005737d-03
		3.59207910759131235356d-03
		1.45620945432529025516d-03
		5.88041240820264096874d-04
		2.46463134818469906812d-04
		7.81794442939557092300d-05
		7.14072491382608190305d-05
		-1.85586374855275456654d-05
		2.59073051863633712884d-05)))

;; kernel tan function on [-pi/4, pi/4], pi/4 ~ 0.7854
;; Input x is assumed to be bounded by ~pi/4 in magnitude.
;; Input y is the tail of x.
;; Input k indicates whether tan (if k = 1) or -1/tan (if k = -1) is returned.
;;
;; Algorithm
;;	1. Since tan(-x) = -tan(x), we need only to consider positive x.
;;	2. if x < 2^-28 (hx<0x3e300000 0), return x with inexact if x!=0.
;;	3. tan(x) is approximated by a odd polynomial of degree 27 on
;;	   [0,0.67434]
;;		  	         3             27
;;	   	tan(x) ~ x + T1*x + ... + T13*x
;;	   where
;;
;; 	        |tan(x)         2     4            26   |     -59.2
;; 	        |----- - (1+T1*x +T2*x +.... +T13*x    )| <= 2
;; 	        |  x 					|
;;
;;	   Note: tan(x+y) = tan(x) + tan'(x)*y
;;		          ~ tan(x) + (1+x*x)*y
;;	   Therefore, for better accuracy in computing tan(x+y), let
;;		     3      2      2       2       2
;;		r = x *(T2+x *(T3+x *(...+x *(T12+x *T13))))
;;	   then
;;		 		    3    2
;;		tan(x+y) = x + (T1*x + (x *(r+y)+y))
;;
;;      4. For x in [0.67434,pi/4],  let y = pi/4 - x, then
;;		tan(x) = tan(pi/4-y) = (1-tan(y))/(1+tan(y))
;;		       = 1 - 2*(tan(y) - (tan(y)^2)/(1+tan(y)))
(declaim (ftype (function (double-float double-float fixnum)
			  double-float)
		kernel-tan))

(defun kernel-tan (x y iy)
  (declare (type (double-float -1d0 1d0) x y)
	   (type (member -1 1) iy)
	   (optimize (speed 3) (safety 0)))
  (let* ((hx (kernel:double-float-high-bits x))
	 (ix (logand hx #x7fffffff))
	 (w 0d0)
	 (z 0d0)
	 (v 0d0)
	 (s 0d0)
	 (r 0d0))
    (declare (double-float w z v s r))
    (when (< ix #x3e300000)
      ;; |x| < 2^-28
      (when (zerop (truncate x))
	(cond ((zerop (logior (logior ix (kernel:double-float-low-bits x))
			      (+ iy 1)))
	       ;; x = 0 and iy = -1 (cot)
	       (return-from kernel-tan (/ (abs x))))
	      ((= iy 1)
	       (return-from kernel-tan x))
	      (t
	       ;; x /= 0 and iy = -1 (cot)
	       ;; Compute -1/(x+y) carefully
	       (let ((a 0d0)
		     (tt 0d0))
		 (setf w (+ x y))
		 (setf z (kernel:make-double-float (kernel:double-float-high-bits w) 0))
		 (setf v (- y (- z x)))
		 (setf a (/ -1 w))
		 (setf tt (kernel:make-double-float (kernel:double-float-high-bits a) 0))
		 (setf s (+ 1 (* tt z)))
		 (return-from kernel-tan (+ tt
					    (* a (+ s (* tt v))))))))))
    (when (>= ix #x3FE59428)
      ;; |x| > .6744
      (when (minusp hx)
	(setf x (- x))
	(setf y (- y)))
      ;; z = pi/4-x
      (setf z (- (kernel:make-double-float #x3FE921FB #x54442D18) x))
      ;; w = pi/4_lo - y
      (setf w (- (kernel:make-double-float #x3C81A626 #x33145C07) y))
      (setf x (+ z w))
      (setf y 0d0))
    (setf z (* x x))
    (setf w (* z z))
    ;; Break x^5*(T[1]+x^2*T[2]+...) into
    ;; x^5(T[1]+x^4*T[3]+...+x^20*T[11]) +
    ;; x^5(x^2*(T[2]+x^4*T[4]+...+x^22*[T12]))
    (setf r (+ (aref ttt 1)
	       (* w
		  (+ (aref ttt 3)
		     (* w
			(+ (aref ttt 5)
			   (* w
			      (+ (aref ttt 7)
				 (* w
				    (+ (aref ttt 9)
				       (* w (aref ttt 11))))))))))))
    (setf v (* z
	       (+ (aref ttt 2)
		  (* w
		     (+ (aref ttt 4)
			(* w
			   (+ (aref ttt 6)
			      (* w
				 (+ (aref ttt 8)
				    (* w
				       (+ (aref ttt 10)
					  (* w (aref ttt 12)))))))))))))
    (setf s (* z x))
    (setf r (+ y (* z (+ (* s (+ r v))
			 y))))
    (incf r (* s (aref ttt 0)))
    (setf w (+ x r))
    (when (>= ix #x3FE59428)
      (let ((v (float iy 1d0)))
	(return-from kernel-tan
	  (* (- 1 (logand 2 (ash hx -30)))
	     (- v
		(* 2
		   (- x (- (/ (* w w)
			      (+ w v))
			   r))))))))
    (when (= iy 1)
      (return-from kernel-tan w))
    ;;
    (let ((a 0d0)
	  (tt 0d0))
      (setf z (kernel:make-double-float (kernel:double-float-high-bits w) 0))
      (setf v (- r (- r x)))		; z + v = r + x
      (setf a (/ -1 w))
      (setf tt (kernel:make-double-float (kernel:double-float-high-bits a) 0))
      (setf s (+ 1 (* tt z)))
      (+ tt
	 (* a
	    (+ s (* tt v)))))))

;; Return sine function of x.
;;
;; kernel function:
;;	__kernel_sin		... sine function on [-pi/4,pi/4]
;;	__kernel_cos		... cose function on [-pi/4,pi/4]
;;	__ieee754_rem_pio2	... argument reduction routine
;;
;; Method.
;;      Let S,C and T denote the sin, cos and tan respectively on 
;;	[-PI/4, +PI/4]. Reduce the argument x to y1+y2 = x-k*pi/2 
;;	in [-pi/4 , +pi/4], and let n = k mod 4.
;;	We have
;;
;;          n        sin(x)      cos(x)        tan(x)
;;     ----------------------------------------------------------
;;	    0	       S	   C		 T
;;	    1	       C	  -S		-1/T
;;	    2	      -S	  -C		 T
;;	    3	      -C	   S		-1/T
;;     ----------------------------------------------------------
;;
;; Special cases:
;;      Let trig be any of sin, cos, or tan.
;;      trig(+-INF)  is NaN, with signals;
;;      trig(NaN)    is that NaN;
;;
;; Accuracy:
;;	TRIG(x) returns trig(x) nearly rounded 
(defun %sin (x)
  (declare (double-float x)
	   (optimize (speed 3)))
  (let ((ix (ldb (byte 31 0) (kernel:double-float-high-bits x))))
    (cond
      ((<= ix #x3fe921fb)
       ;; |x| < pi/4, approx
       (kernel-sin x 0d0 0))
      ((>= ix #x7ff00000)
       ;; sin(Inf or NaN) is NaN
       (- x x))
      (t
       ;; Argument reduction needed
       (multiple-value-bind (n y0 y1)
	   (kernel::%ieee754-rem-pi/2 x)
	 (case (logand n 3)
	   (0
	    (kernel-sin y0 y1 1))
	   (1
	    (kernel-cos y0 y1))
	   (2
	    (- (kernel-sin y0 y1 1)))
	   (3
	    (- (kernel-cos y0 y1)))))))))

(defun %cos (x)
  (declare (double-float x)
	   (optimize (speed 3)))
  (let ((ix (ldb (byte 31 0) (kernel:double-float-high-bits x))))
    (cond
      ((< ix #x3fe921fb)
       ;;|x| < pi/4, approx
       (kernel-cos x 0d0))
      ((>= ix #x7ff00000)
       ;; cos(Inf or NaN) is NaN
       (- x x))
      (t
       ;; Argument reduction needed
       (multiple-value-bind (n y0 y1)
	   (kernel::%ieee754-rem-pi/2 x)
	 (ecase (logand n 3)
	   (0
	    (kernel-cos y0 y1))
	   (1
	    (- (kernel-sin y0 y1 1)))
	   (2
	    (- (kernel-cos y0 y1)))
	   (3
	    (kernel-sin y0 y1 1))))))))

(defun %tan (x)
  (declare (double-float x)
	   (optimize (speed 3)))
  (let ((ix (logand #x7fffffff (kernel:double-float-high-bits x))))
    (cond ((<= ix #x3fe921fb)
	   (kernel-tan x 0d0 1))
	  ((>= ix #x7ff00000)
	   (- x x))
	  (t
	   (multiple-value-bind (n y0 y1)
	       (kernel::%ieee754-rem-pi/2 x)
	     (let ((flag (- 1 (ash (logand n 1) 1))))
	       ;; flag = 1 if n even, -1 if n odd
	       (kernel-tan y0 y1 flag)))))))

(declaim (ext:end-block))