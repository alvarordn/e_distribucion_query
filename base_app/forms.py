from django.forms import ModelForm
from .models import EndesaCredentials

class EndesaCredentialsForm(ModelForm):    
    class Meta:
        model = EndesaCredentials
        fields = ['username', 'password', 'date']

