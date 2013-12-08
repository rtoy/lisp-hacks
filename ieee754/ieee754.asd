;; -*- Mode: lisp -*-

(asdf:defsystem ieee754
  :author "Raymond Toy"
  :components
  ((:file "package")
   (:file "k-rem-pi")
   (:file "e-rem-pi" :depends-on ("k-rem-pi"))))

(defmethod perform ((op asdf:test-op) (c (eql (asdf:find-system :ieee754))))
  (asdf:oos 'asdf:test-op :ieee754-tests))

(asdf:defsystem ieee754-tests
  :depends-on (ieee754)
  :in-order-to ((compile-op (load-op :rt))
		(test-op (load-op :rt :oct)))
  :components
  ((:file "test")))

(defmethod perform ((op asdf:test-op) (c (eql (asdf:find-system :ieee754-tests))))
  (funcall (intern "DO-TESTS" (find-package "RT"))))

