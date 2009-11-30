;; test the smash alg i put together

(defun smash (region results-list)
	(let* ((local-results (split region))
				 (n (length local-results)))
		(cond ((= n 0) (format t "error no regions"))
					((= n 1) (setf results-list (append results-list local-results)))
					((> n 1) (dotimes (i n)
										 (setf results-list (append results-list (list (smash (nth i local-results) results-list)))))
					 results-list))))


(defun split (region)
	(let ((prob (get-random-uniform)))
		(if (> 0.7 prob) (list 'A) (list 'B 'C))))
					


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

				