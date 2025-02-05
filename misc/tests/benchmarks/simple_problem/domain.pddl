(define (domain finite-domain)
  (:requirements :strips)
  (:predicates (value ?v - object ?val - object))

  (:action o1
    :parameters ()
    :precondition (value a 0)
    :effect (and (not (value a 0)) (value a 1)))

  (:action o2
    :parameters ()
    :precondition (value b 0)
    :effect (and (not (value b 0)) (not (value a 1)) (value a 0) (value b 1))))
