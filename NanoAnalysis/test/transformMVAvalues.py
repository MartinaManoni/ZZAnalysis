from decimal import Decimal, getcontext
import math

# Set precision (increase if needed)
getcontext().prec = 10

# Input values
values = [
    Decimal('1.633973689084034'),
    Decimal('1.5499076306249353'),
    Decimal('2.0629564440753247'),
    Decimal('0.3685228146685872'),
    Decimal('0.2662407818935475'),
    Decimal('-0.5444837363886459')
]

def transform(x: Decimal) -> Decimal:
    exp_part = Decimal(math.exp(-2 * float(x)))
    return Decimal(2) / (Decimal(1) + exp_part) - Decimal(1)

# Transform and print results
for val in values:
    result = transform(val)
    print(f"Input: {val} -> Transformed: {result}")