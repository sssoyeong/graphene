from django.shortcuts import render
from django.views.generic import TemplateView
# from forms import AtomForm
from . import forms
from . import predict
import numpy as np

class IndexView(TemplateView):
    template_name = "ml/index.html"
    form_class = forms.AtomForm

    def get_context_data(self,  **kwargs):
        context = super(IndexView, self).get_context_data(**kwargs)
        context['atomRange'] = range(0, 36)
        context['atomNum'] = range(0,2)
        context['form'] = self.form_class

        return context

    def post(self, request, *args, **kwargs):
        form = self.form_class(request.POST)
        context = {}
        context['atomRange'] = range(0, 36)
        context['atomNum'] = range(0,2)
        context['form'] = self.form_class

        if form.is_valid():
            atom_category = form.cleaned_data['atomCategory']
            atom_list = form.cleaned_data['atomList']
            broad  = float(form.cleaned_data['broad'])
            erange = float(form.cleaned_data['erange'])
            predicted_data = predict.predict(atom_category,atom_list, broad, erange)

            #result_y = []
            #start = -10.0
            #while start <= 10.0:
            #    result_y.append(start)
            #    start = start + 0.02
            result_y = list( np.linspace(-erange, erange, 1001) )

            final_result = []
            for i in range(len(result_y)):
                temp = []
                temp.append(result_y[i])
                temp.append(predicted_data[0][i])
                temp.append(predicted_data[1][i])
                final_result.append(temp)
            context['result'] = final_result
        return render(request, self.template_name, context)
