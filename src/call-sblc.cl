(cl:defpackage "CALL-EMULATOR" (:use "CL" "SB-ALIEN" "SB-C-CALL"))
(cl:in-package "CALL-EMULATOR")

(defparameter *nmodel-points* 10)
(defparameter *nemu-points* 40)

;;; Define the lisp function interface to the C routine
;;; not sure about the declaim call
(declaim (inline callEmulator))
(define-alien-routine callEmulator
		(void)
	(xmodel-in (* (array double *nmodel-points*)))
	(nparams-in (* int))
	(training-in (* (array double *nmodel-points*)))
	(nmodelpts (* int))
	(nthetas-in (* int))
	(final-emulated-x (* (array double *nemu-points*)))
	(nemupts-in (* int))
	(final-emulated-y (* (array double *nemu-points*)))
	(final-emulated-var (* (array double *nemu-points*)))
	(range-min  (* double))
	(range-max (* double)))

;;;  a function which sets up the values we want to pass 
;;; and then actually calls the emulator
(defun call-emulator (xmodel training)
	...