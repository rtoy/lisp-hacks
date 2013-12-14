;; -*- Mode: lisp -*-

(asdf:defsystem fdlibm
  :author "Raymond Toy"
  :components
  ((:file "fdlibm")
   (:file "trig" :depends-on ("fdlibm"))))

(asdf:defsystem fdlibm-tests
  :author "Raymond Toy"
  :depends-on (fdlibm)
  :in-order-to ((compile-op (load-op :rt))
		(test-op (load-op :rt)))
  :components
  ((:file "fdlibm-tests")))

(defmethod perform ((op asdf:test-op) (c (eql (asdf:find-system :fdlibm-tests))))
  (funcall (intern "DO-TESTS" (find-package "RT"))))

