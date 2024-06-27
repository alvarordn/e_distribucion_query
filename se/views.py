from rest_framework.decorators import api_view
from rest_framework.response import Response
from .serializers import StateEstimatorSerializer
from base_app.models import StateEstimator
from .lib import grid 


@api_view(['GET'])
def getProjects(request):
    Projects = StateEstimator.objects.all()
    serializer = StateEstimatorSerializer(Projects, many=True)
    return Response(serializer.data)

@api_view(['GET'])
def getProject(request, pk):
    Projects = StateEstimator.objects.get(id=pk)
    serializer = StateEstimatorSerializer(Projects, many=False)
    return Response(serializer.data)

@api_view(['POST'])
def solveProblem(request, pk):
    project = StateEstimator.objects.get(id=pk)
   
    nodes = eval(project.nodes)
    lines = eval(project.lines)
    meas = eval(project.meas)
    net = grid(nodes, lines, meas)
    x = [0, 0, 1, 1, 1]
    res, sol = net.state_estimation(x)
    print(sol[-1])
    
    serializer = StateEstimatorSerializer(project, many=False)
    return Response(serializer.data)