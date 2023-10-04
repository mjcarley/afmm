(TeX-add-style-hook
 "notes"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "amsmath"
    "times")
   (TeX-add-symbols
    '("gnbarijk" 4)
    '("gnijk" 4)
    "D"
    "E"
    "I"
    "bA")
   (LaTeX-add-labels
    "equ:basic"
    "equ:gfunc"
    "equ:recursion"
    "equ:gfunc:dx"))
 :latex)

