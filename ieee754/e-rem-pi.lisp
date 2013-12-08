(in-package "IEEE754")

(defconstant 2/pi-bits
  (make-array 66 :element-type '(unsigned-byte 32)
	      :initial-contents
	      '(#xA2F983 #x6E4E44 #x1529FC #x2757D1 #xF534DD #xC0DB62 
		#x95993C #x439041 #xFE5163 #xABDEBB #xC561B7 #x246E3A 
		#x424DD2 #xE00649 #x2EEA09 #xD1921C #xFE1DEB #x1CB129 
		#xA73EE8 #x8235F5 #x2EBB44 #x84E99C #x7026B4 #x5F7E41 
		#x3991D6 #x398353 #x39F49C #x845F8B #xBDF928 #x3B1FF8 
		#x97FFDE #x05980F #xEF2F11 #x8B5A0A #x6D1F6D #x367ECF 
		#x27CB09 #xB74F46 #x3F669E #x5FEA2D #x7527BA #xC7EBE5 
		#xF17B3D #x0739F7 #x8A5292 #xEA6BFB #x5FB11F #x8D5D08 
		#x560330 #x46FC7B #x6BABF0 #xCFBC20 #x9AF436 #x1DA9E3 
		#x91615E #xE61B08 #x659985 #x5F14A0 #x68408D #xFFD880 
		#x4D7327 #x310606 #x1556CA #x73A8C9 #x60E27B #xC08C6B))
  "Table of constants for 2/pi, 396 hex digits")

(defconstant npio2-hw
  (make-array 32 :element-type '(unsigned-byte 32)
	      :initial-contents
	      '(#x3FF921FB #x400921FB #x4012D97C #x401921FB #x401F6A7A #x4022D97C
		#x4025FDBB #x402921FB #x402C463A #x402F6A7A #x4031475C #x4032D97C
		#x40346B9C #x4035FDBB #x40378FDB #x403921FB #x403AB41B #x403C463A
		#x403DD85A #x403F6A7A #x40407E4C #x4041475C #x4042106C #x4042D97C
		#x4043A28C #x40446B9C #x404534AC #x4045FDBB #x4046C6CB #x40478FDB
		#x404858EB #x404921FB)))
  

(defconstant two24 (scale-float 1d0 24))
(defconstant 2/pi (/ 2 pi))
(defconstant pi/2-1  (kernel:make-double-float #x3FF921FB #x54400000))
(defconstant pi/2-1t (kernel:make-double-float #x3DD0B461 #x1A626331))
(defconstant pi/2-2  (kernel:make-double-float #x3DD0B461 #x1A600000))
(defconstant pi/2-2t (kernel:make-double-float #x3BA3198A #x2E037073))
(defconstant pi/2-3  (kernel:make-double-float #x3BA3198A #x2E000000))
(defconstant pi/2-3t (kernel:make-double-float #x397B839A #x252049C1))

(defun ieee754-rem-pi/2 (x yy)
  (declare (double-float x)
	   (type (simple-array double-float (*)) yy)
	   (optimize (speed 3)))
  (let* ((result-n 0)
	 (result-y0 0d0)
	 (result-y1 0d0)
	 (z 0d0)
	 (hx (kernel:double-float-high-bits x))
	 (ix (logand hx #x7fffffff)))
    (declare (double-float z result-y0 result-y1))
    #+nil
    (format t "ix = ~A~%" ix)
    (cond ((<= ix #x3fe921fb)
	   ;; |x| <= pi/4, no need for reduction
	   (setf result-n 0
		 result-y0 x
		 result-y1 0d0))
	  ((< ix #x4002d97c)
	   ;; |x| < 3pi/4. Special case with n = +/- 1
	   (cond ((plusp hx)
		  (setf z (- x pi/2-1))
		  (cond ((/= ix #x3ff921fb)
			 ;; 33 + 53 bits of pi is enough
			 (let ((y0 (- z pi/2-1t)))
			   (setf result-n 1
				 result-y0 y0
				 result-y1 (- (- z y0) pi/2-1t))))
			(t
			 ;; near pi/2. Use 33+33+53 bits of pi
			 (format t "near pi/2~%")
			 (decf z pi/2-2)
			 (let ((y0 (- z pi/2-2t)))
			   (setf result-n 1
				 result-y0 y0
				 result-y1 (- (- z y0) pi/2-2t))))))
		 (t
		  ;; Negative x
		  (setf z (+ x pi/2-1))
		  (cond ((/= ix #x3ff921fb)
			 ;; 33 + 53 bits of pi is enough
			 (let ((y0 (+ z pi/2-1t)))
			   (setf result-n -1
				 result-y0 y0
				 result-y1 (+ (- z y0) pi/2-1t))))
			(t
			 ;; near pi/2. Use 33+33+53 bits of pi
			 (incf z pi/2-2)
			 (let ((y0 (+ z pi/2-2t)))
			   (setf result-n -1
				 result-y0 y0
				 result-y1 (+ (- z y0) pi/2-2t))))))))
	  ((<= ix #x413921fb)
	   ;; |x| <= 2^19*pi/2, medium size
	   (let* ((tt (abs x))
		  (n (truncate (+ (* tt 2/pi) 0.5d0)))
		  (fn (float n 1d0))
		  (r (- tt (* fn pi/2-1)))
		  (w (* fn pi/2-1t))
		  (y0 0d0)
		  (y1 0d0))
	     (declare (type (double-float 0d0 1048576d0) tt)
		      (double-float r w y0 y1))
	     ;;#+nil
	     (format t "n, fn = ~S ~S~%" n fn)
	     ;; First round good to 85 bit
	     (cond ((and (< fn 32)
			 (/= ix (aref npio2-hw (- n 1))))
		    ;;#+nil
		    (format t "first round ~%")
		    (setf y0 (- r w)))
		   (t
		    (setf y0 (- r w))
		    ;; j is the IEEE exponent of x. The logand part is
		    ;; the IEEE exponent of y0.  Can we do this
		    ;; faster?
		    (let* ((j (ash ix -20))
			   (i (- j (logand (ash (kernel:double-float-high-bits y0) -20)
					   #x7ff))))
		      ;;#+nil
		      (format t "i = ~S~%" i)
		      (when (> i 16)
			;; 2nd iteration, good to 118
			#+nil
			(format t "second round~%")
			(setf tt r)
			(setf w (* fn pi/2-2))
			(setf r (- tt w))
			(setf w (- (* fn pi/2-2t)
				   (- (- tt r)
				      w)))
			(setf y0 (- r w))
			(let ((i (- j (logand (ash (kernel:double-float-high-bits y0) -20)
					      #x7ff))))
			  #+nil
			  (format t "i = ~S~%" i)
			  (when (> i 49)
			    ;; 3rd iteration needed. 151 bits
			    #+nil
			    (format t "third round~%") 
			    (setf tt r)
			    (setf w (* fn pi/2-3))
			    (setf r (- tt w))
			    (setf w (- (* fn pi/2-3t)
				       (- (- tt r)
					  w)))
			    (setf y0 (- r w))))))))
	     (setf y1 (- (- r y0) w))
	     (cond ((minusp hx)
		    (setf result-n (- n)
			  result-y0 (- y0)
			  result-y1 (- y1)))
		   (t
		    (setf result-n n
			  result-y0 y0
			  result-y1 y1)))))
	  ((>= ix #x7ff00000)
	   ;; x is inf or NaN
	   (let ((y (- x x)))
	   (setf result-n 0
		 result-y0 y
		 result-y1 y)))
	  (t
	   ;; All other large values
	   #+nil
	   (format t "Big value ~S~%" x)
	   (let* ((ilogx (ext:truly-the (integer 0 1023) (kernel::logb (abs x))))
		  (z (scale-float (abs x) (- 23 ilogx)))
		  (e0 (- ilogx 23))
		  (y (make-array 2 :element-type 'double-float))
		  (tx (make-array 3 :element-type 'double-float))
		  (nx 3))
	     (declare (type (integer 0 3) nx))
	     #+nil
	     (format t "z = ~S~%" z)
	     (dotimes (i 2)
	       (setf (aref tx i) (ftruncate (the (double-float 0d0 1048576d0) z)))
	       (setf z (* (- z (aref tx i)) (scale-float 1d0 24))))
	     (setf (aref tx 2) z)
	     #+nil
	     (format t "tx = ~S~%" tx)
	     (loop while (zerop (aref tx (- nx 1)))
		   do (decf nx))
	     #+nil
	     (progn
	       (format t "nx = ~S~%" nx)
	       (format t "e0 = ~S~%" e0))
	     (let ((n (kernel-rem-pi/2 tx y e0 nx 2 2/pi-bits)))
	       (cond ((minusp hx)
		      (setf result-n (- n)
			    result-y0 (- (aref y 0))
			    result-y1 (- (aref y 1))))
		     (t
		      (setf result-n n
			    result-y0 (aref y 0)
			    result-y1 (aref y 1))))))))
    (setf (aref yy 0) result-y0)
    (setf (aref yy 1) result-y1)
    (values result-n)))
