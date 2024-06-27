from django.shortcuts import render, redirect
from .forms import EndesaCredentialsForm
from .utils import procesar_credenciales

# Create your views here.
def base_view(request):
    return render(request, 'base_app/base_view.html')


def endesa(request):
    resultado = None 
    if request.method == 'POST':
        ECform = EndesaCredentialsForm({'username': request.POST['username'],
                              'password': request.POST['password'],
                              'date': request.POST['date']})      
        if ECform.is_valid():
            username = ECform.cleaned_data['username']
            password = ECform.cleaned_data['password']
            date = ECform.cleaned_data['date']
            ECform.save()
            resultado = procesar_credenciales(username, password, date)
    else:
        ECform = EndesaCredentialsForm()
    return render(request, 'base_app/query.html', {'form': ECform, 'resultado': resultado})
