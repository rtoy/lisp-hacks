(in-package "IEEE754")

#+nil
(progn
(declaim (inline rem-pi/2))
(alien:def-alien-routine ("rem_pio2" rem-pi/2) c-call:int
  (x double-float)
  (y0 double-float :out)
  (y1 double-float :out))
)

(defun compare (x)
  (let* ((y (make-array 2 :element-type 'double-float))
	 (n (ieee754-rem-pi/2 x y)))
    (multiple-value-bind (true-n true-y0 true-y1)
	(kernel::%ieee754-rem-pi/2 x)
      (unless (and (= n true-n)
		   (eql (aref y 0) true-y0)
		   (eql (aref y 1) true-y1))
	(list x
	      (list n (aref y 0) (aref y 1))
	      (list true-n true-y0 true-y1))))))

(defun test-mult-pi/4 (n &optional (count 1000))
  ;; Test multiples of random multiples of pi up to n*pi
  (let ((pi/4 (/ pi 4))
	(fail-count 0))
    (dotimes (k count)
      (let* ((x (* (random n) pi/4))
	     (result (compare x)))
	(when result
	  (incf fail-count)
	  (format t "~S~%" result))))
    (format t "~D failures out of ~D tests~%" fail-count count)))

(defun test-time (x n)
  (declare (double-float x)
	   (fixnum n))
  (flet ((test-lisp ()
	   (let ((nn 0)
		 (yy0 0d0)
		 (yy1 0d0)
		 (yy (make-array 2 :element-type 'double-float)))
	     (declare (fixnum n)
		      (double-float yy0 yy1))
	     (dotimes (k n)
	       (declare (fixnum k))
	       (multiple-value-bind (n)
		   (ieee754-rem-pi/2 x yy)
		 (setf nn n)
		 (incf yy0 (aref yy 0))
		 (incf yy1 (aref yy 1))))
	     (values nn yy0 yy1)))
	 (test-ref ()
	   (let ((nn 0)
		 (yy0 0d0)
		 (yy1 0d0))
	     (declare (fixnum n)
		      (double-float yy0 yy1))
	     (dotimes (k n)
	       (declare (fixnum k))
	       (multiple-value-bind (n y0 y1)
		   (kernel::%ieee754-rem-pi/2 x)
		 (setf nn n)
		 (incf yy0 y0)
		 (incf yy1 y1)))
	     (values nn yy0 yy1))))
    (multiple-value-bind (n y0 y1)
	(time (test-lisp))
      (format t "Actual:   ~S ~S ~S~%" n y0 y1))
    (multiple-value-bind (n y0 y1)
	(time (test-ref))
      (format t "Expected: ~S ~S ~S~%" n y0 y1))))

(rt:deftest small.1
    (compare .1d0)
  nil)

(rt:deftest nsmall.1
    (compare -.1d0)
  nil)

(rt:deftest small.2
    (compare .5d0)
  nil)

(rt:deftest nsmall.2
    (compare -.5d0)
  nil)

(rt:deftest pi/4.1
    (compare (/ pi 4))
  nil)

(rt:deftest npi/4.1
    (compare (/ pi -4))
  nil)

(rt:deftest pi.1
    (compare pi)
  nil)

;; Tests the case |x| < 3*pi/4; Extra bits needed near pi/2
(rt:deftest pi/2.1
    (compare (/ pi 2))
  nil)

;; Tests the case |x| < 3*pi/4;  33+53 bits is enough
(rt:deftest special.1
    (compare 1.5d0)
  nil)

(rt:deftest special.2
    (compare 2d0)
  nil)

(rt:deftest special.3
    (compare (* 3/4 pi))
  nil)

;; Medium size 3/4*pi < x < 2^19*pi/2
(rt:deftest medium.1
    (compare pi)
  nil)

;; All of the cases in npio2-hw, first round case
(rt:deftest medium-first-round.1
    (let (results)
      (loop for x across ieee754::npio2-hw
	    do
	       (push (compare (* 1.01d0 (kernel:make-double-float x 0)))
		     results))
      (remove nil results))
  nil)

;; All of the cases in npio2-hw, not first round case. Second round
;; case is covered.
(rt:deftest medium-second-round.1
    (let (results)
      (loop for x across ieee754::npio2-hw
	    do
	       (push (compare (kernel:make-double-float x 0))
		     results))
      (remove nil results))
  nil)

;; All of the cases in npio2-hw, not first round case. Third round
;; case is covered.
(rt:deftest medium-third-round.1
    (compare (* 5 pi))
  nil)


(rt:deftest inf.1
    (ext:with-float-traps-masked (:invalid)
      (compare ext:double-float-positive-infinity))
  nil)

(rt:deftest nan.1
    (ext:with-float-traps-masked (:invalid)
      (compare (kernel:make-double-float #x7ff00000 1)))
  nil)

;; Tests for values bigger than 2^19*pi/2.  These call k-rem-pi.
(rt:deftest basic.big.1
    (compare (scale-float 1d0 120))
  nil)

(rt:deftest basic.big.2
    (compare (scale-float 1d0 1023))
  nil)

(rt:deftest basic.big.3
    (compare most-positive-double-float)
  nil)
(rt:deftest basic.big.4
    (compare most-negative-double-float)
  nil)
