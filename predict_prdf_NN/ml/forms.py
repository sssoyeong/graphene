from django import forms
ATOM_CHOICES = (
    ('B', 'B'),
    ('N', 'N'),
)

ENERGY_CHOICES = (
    (10, 10),
    (5, 5),
    (3, 3),
)

BROAD_CHOICES = (
    (0.05, 0.05),
    (0.10, 0.10),
)

class AtomForm(forms.Form):
    atomList = forms.CharField()
    atomCategory = forms.ChoiceField(choices=ATOM_CHOICES)
    broad  = forms.ChoiceField(choices=BROAD_CHOICES)
    erange = forms.ChoiceField(choices=ENERGY_CHOICES)
