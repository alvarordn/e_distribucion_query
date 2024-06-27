from rest_framework import serializers
from base_app.models import StateEstimator

class StateEstimatorSerializer(serializers.ModelSerializer):
    class Meta:
        model = StateEstimator
        fields = '__all__'