(in-package "CL-USER")

(defconstant init-jk
  (make-array 4 :element-type '(unsigned-byte 8)
	      :initial-contents '(2 3 4 6)))

(defconstant pio2
  (make-array 8 :element-type 'double-float
	      :initial-contents (list
				 (kernel:make-double-float #x3FF921FB #x40000000)
				 (kernel:make-double-float #x3E74442D #x00000000)
				 (kernel:make-double-float #x3CF84698 #x80000000)
				 (kernel:make-double-float #x3B78CC51 #x60000000)
				 (kernel:make-double-float #x39F01B83 #x80000000)
				 (kernel:make-double-float #x387A2520 #x40000000)
				 (kernel:make-double-float #x36E38222 #x80000000)
				 (kernel:make-double-float #x3569F31D #x00000000))))

(defun kernel-rem-pi/2 (x y e0 nx prec ipio2)
  (declare (type (simple-array double-float (*)) x y)
	   (type (integer 0 3) prec)
	   (type (simple-array (unsigned-byte 32) (*)) ipio2))
  (let* ((jk (aref init-jk prec))
	 (jp jk)
	 (jx (- nx 1))
	 (jv (truncate (- e0 3) 24))
	 (q0 0)
	 (iq (make-array 20 :element-type '(signed-byte 32)))
	 (f (make-array 20 :element-type 'double-float))
	 (fq (make-array 20 :element-type 'double-float))
	 (q (make-array 20 :element-type 'double-float)))
    (when (minusp jv)
      (setf jv 0))
    (setf q0 (- e0 (* 24 (+ jv 1))))
    #+nil
    (format t "q0 = ~A~%" q0)
    ;; set up f[0] to f[jx+jk] where f[jx+jk] = ipio2[jv+jk]
    (let ((j (- jv jx))
	  (m (+ jx jk)))
      #+nil
      (format t "j = ~S, m = ~S~%" j m)
      (loop for i from 0 upto m
	    do
	       (progn
		 (setf (aref f i) (if (minusp j) 0d0 (float (aref ipio2 j) 1d0)))
		 (incf j))))
    #+nil
    (format t "f = ~S~%" f)
    ;; Compute q[0],...,q[jk]
    (loop for i from 0 upto jk do
      (let ((fw 0d0))
	(loop for j from 0 upto jx do
	  (progn
	    (incf fw (* (aref x j) (aref f (- (+ jx i) j))))
	    (setf (aref q i) fw)))))
    #+nil
    (format t "q = ~S~%" q)
    (let ((jz jk)
	  (n 0)
	  (ih 0))
      (tagbody
	recompute
	 ;; distill q[] into iq[] reversingly
	 (let ((i 0)
	       (z (aref q jz))
	       (fw 0d0))
	   (loop for j from jz above 0
		 do
		    (let ((fw (ftruncate (* z (scale-float 1d0 -24)))))
		      (setf (aref iq i) (truncate (- z (* fw (scale-float 1d0 24)))))
		      #+nil
		      (format t "i ~D: z = ~S fw = ~S iq = ~S~%"
			      i z fw (aref iq i))
		      (setf z (+ fw (aref q (- j 1))))
		      (incf i)))
	   #+nil
	   (format t "init iq = ~S~%" iq)
	   ;; compute n
	   (setf z (scale-float z q0))
	   (setf z (- z (* 8 (floor (* z 0.125d0)))))
	   (setf n (truncate z))
	   (setf z (- z n))
	   (setf ih 0)
	   #+nil
	   (format t "q0 = ~S, z = ~S~%" q0 z)
	   (cond ((plusp q0)
		  ;; need iq[jz-1] to determine n
		  (setf i (ash (aref iq (- jz 1)) (- q0 24)))
		  (incf n i)
		  (decf (aref iq (- jz 1)) (ash i (- 24 q0)))
		  (setf ih (ash (aref iq (- jz 1)) (- q0 23))))
		 ((zerop q0)
		  (setf ih (ash (aref iq (- jz 1)) -23)))
		 ((>= z 0.5)
		  (setf ih 2)))
	   #+nil
	   (format t "iq = ~S~%" iq)
	   (when (plusp ih)
	     ;; q > 0.5
	     #+nil
	     (format t "q > 0.5~%")
	     (incf n 1)
	     (let ((carry 0)
		   (j 0))
	       (loop for i from 0 below jz
		     do
			(progn
			  (setf j (aref iq i))
			  (cond ((zerop carry)
				 (when (/= j 0)
				   (setf carry 1)
				   (setf (aref iq i) (- #x1000000 j))))
				(t
				 (setf (aref iq i) (- #xffffff j))))))
	       #+nil
	       (progn
		 (format t "iq = ~S~%" iq)
		 (format t "q = ~S~%" q))
	       (when (plusp q0)
		 ;; rare case: chance is 1 in 12
		 #+nil
		 (format t "rare case of q0 > 0~%")
		 (case q0
		   (1
		    (setf (aref iq (- jz 1))
			  (logand (aref iq (- jz 1)) #x7fffff)))
		   (2
		    (setf (aref iq (- jz 1))
			  (logand (aref iq (- jz 1)) #x3fffff)))))
	       (when (= ih 2)
		 (setf z (- 1 z))
		 (when (/= carry 0)
		   (decf z (scale-float 1d0 q0))))))
	   ;; Check if recomputation is needed
	   #+nil
	   (format t "recomp? z = ~S~%" z)
	   (when (zerop z)
	     #+nil
	     (format t "Checking if recomp needed~%")
	     (let ((j 0))
	       (loop for i from (- jz 1) downto jk
		     do (setf j (logior j (aref iq i))))
	       (when (zerop j)
		 ;; need recomputation
		 ;; k = no of terms needed
		 #+nil
		 (format t "iq = ~S~%" iq)
		 (let ((k
			(loop for k from 1 while (zerop (aref iq (- jk k)))
			      finally (return k))))
		   #+nil
		   (format t "k = ~S~%" k)
		   (loop for i from (+ jz 1) upto (+ jz k)
			 do
			    (progn
			      ;; add q[jz + 1] to q[jz + k]
			      (setf (aref f (+ jx i)) (float (aref ipio2 (+ jv i)) 1d0))
			      (loop for j from 0 upto jx with fw = 0d0
				    do
				       (incf fw (* (aref x j) (aref f (- (+ jx i) j))))
				    finally
				       (setf (aref q i) fw))))
		   (incf jz k)
		   #+nil
		   (format t "recompute needed! jz = ~D~%" jz)
		   (go recompute)))))

	   ;; chop off zero terms
	   #+nil
	   (format t "chop off zero: z = ~S~%" z)
	   (cond ((zerop z)
		  (decf jz)
		  (decf q0 24)
		  (loop while (zerop (aref iq jz))
			do
			   (progn
			     (decf jz)
			     (decf q0 24))))
		 (t
		  ;; break z into 24-bit if necessary
		  (setf z (scale-float z (- q0)))
		  (cond ((>= z (scale-float 1d0 24))
			 (setf fw (ftruncate (* z (scale-float 1d0 -24))))
			 (setf (aref iq jz) (truncate (- z
							 (* fw (scale-float 1d0 24)))))
			 (incf jz)
			 (incf q0 24)
			 (setf (aref iq jz) (truncate fw)))
			(t
			 (setf (aref iq jz) (truncate z))))))
	   ;; convert integer "bit" chunk to floating-point value
	   (setf fw (scale-float 1d0 q0))
	   (loop for i from jz downto 0
		 do
		    (progn
		      (setf (aref q i) (* fw (aref iq i)))
		      (setf fw (* fw (scale-float 1d0 -24)))))
	   ;; Compute PIo2[0,....,jp]*q[jz,...,0]
	   (loop for i from jz downto 0
		 do
		    (loop with fw = 0d0
			  for k from 0
			  while (and (<= k jp) (<= k (- jz i)))
			  do
			     (incf fw (* (aref pio2 k) (aref q (+ i k))))
			  finally
			     (setf (aref fq (- jz i)) fw)))

	   ;; Compress fq[] in to y[]
	   (ecase prec
	     (0
	      (setf fw 0d0)
	      (loop for i from jz downto 0
		    do (incf fw (aref fq i)))
	      (setf (aref y 0) (if (zerop ih) fw (- fw))))
	     ((1 2)
	      (setf fw 0d0)
	      (loop for i from jz downto 0
		    do (incf fw (aref fq i)))
	      (setf (aref y 0) (if (zerop ih) fw (- fw)))
	      (setf fw (- (aref fq 0) fw))
	      (loop for i from 1 upto jz
		    do (incf fw (aref fq i)))
	      (setf (aref y 1) (if (zerop ih) fw (- fw)))))))
      (logand n 7))))

    
	   
