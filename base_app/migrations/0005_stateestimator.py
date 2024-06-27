# Generated by Django 5.0.6 on 2024-06-11 10:27

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('base_app', '0004_powerflow_lines_powerflow_loads_powerflow_nodes'),
    ]

    operations = [
        migrations.CreateModel(
            name='StateEstimator',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reference', models.CharField(max_length=200)),
                ('name', models.CharField(max_length=200)),
                ('email', models.CharField(max_length=200)),
                ('rol', models.CharField(choices=[('1', 'Student'), ('2', 'Professor'), ('3', 'Other')], default=None, max_length=200)),
                ('nodes', models.TextField(default=' ', max_length=2000)),
                ('lines', models.TextField(default=' ', max_length=2000)),
                ('meas', models.TextField(default=' ', max_length=2000)),
            ],
        ),
    ]