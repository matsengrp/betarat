"""Implementation of Hyp2F1 implemented in C using a finitely terminating, recursively factored form of the
function, as seen in Grahm, Knuth, Patashnik's Concrete Mathematics, section 5.5"""

# This makes the function we want accessible via:
#     from betrat import rf_hyp2f1
#     rf_hyp2f1.hyp2f1(...)
from rf_hyp2f1 import hyp2f1
