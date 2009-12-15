
;; returns a float in the range a-> b
;; this is not portable, i have "guessed" that the max int is 2^32
;; this may not be true even on  32 bit systems and certainly not on
;; 64bit, what if random -> long int? 
;; i don't know how to deal with this yet, but testing seems to be required
(defun get-random-float-in-range (A B) 
	(let ((flt-max (expt 2 32)) (my-random 0))
		(setf my-random (float (/ (random flt-max) flt-max)))
		;;debug (format t "~a%" my-random)
		(setf my-random (+ (* my-random (- B A)) A))
		my-random))

(defun get-random-uniform () 
	(get-random-float-in-range 0 1))

(defun split (region depth)
	(let ((p (get-random-uniform)))
		(if (< 0.2 p) (list depth) (list depth depth))))
	
(defun smash (region region-list depth)
	(let* ((local-list (split region depth))
				 (n (length local-list)))
		(cond ((= n 0) 'ERROR)
					((= n 1) (setf region-list (append region-list local-list)))
					((> n 1) (dotimes (i n)
										 (setf region-list 
													 (append region-list 
																	 (smash (nth i local-list ) region-list (incf depth)))))))))
