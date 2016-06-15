# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 18:07:35 2016

@author: agniv
"""
import numpy
import math
import transformations as tr

def getX(x):
    return numpy.matrix(
    [[1,0,0],
     [0,math.cos(x),-math.sin(x)],
     [0,math.sin(x),math.cos(x)]])
     

def getY(y):
    return numpy.matrix(
    [[math.cos(y),0,math.sin(y)],
     [0,1,0],
     [-math.sin(y),0,math.cos(y)]])    
    
def getZ(z):
    return numpy.matrix(
    [[math.cos(z),-math.sin(z),0],
     [math.sin(z),math.cos(z),0],
     [0,0,1]])    
    

def rotationMatrix(roll, pitch, yaw, seq):
    """Form a 3x3 rotation matrix from the given Euler angles.
    
    Keyword arguments:
    roll -- rotation along X-axis
    pitch -- rotation along Y-axis
    yaw -- rotation along Z-axis
    seq -- three digit number denoting the order of rotation, in X-Y-Z sequence 
            needs to be a three digit number, with each number denoting the
            relative rank of that axis
            Example:'312' denotes: rotation along Y, followed by Z, followed by X
                    '132' denotes: rotation along X, followed by Z, followed by Y
    """
    seq0 = seq%10
    seq1 = ((seq - seq0)/10)%10
    seq2 = ((((seq - seq0)/10) - seq1)/10)%10
    
    if ((seq0==seq1==seq2) or ((seq0+seq1+seq2)!=6) or ((seq0*seq1*seq2)!=6) or ((seq0-1)*(seq1-1)*(seq2-1)!=0)):
        raise ValueError('The sequence string has been formatted wrong!')
    else:
        rotX = getX(roll)
        rotY = getY(pitch)
        rotZ = getZ(yaw)
        A = rotX
        B = rotY
        C = rotZ
        
        if(seq1==1):
            A = rotY
            if(seq2==2):
                B = rotX
                C = rotZ
            else:
                B = rotZ
                C = rotX       
        elif(seq0==1):
            A = rotZ
            if(seq2==2):
                B = rotX
                C = rotY
            else:
                B = rotY
                C = rotX
        elif(seq2==1):
            if(seq1==3):
                B = rotZ
                C = rotY
        
        return numpy.dot(numpy.dot(C,B),A)


def eulerLogarithm(rotationMatrix):
    """Form the 3x1 Lie algebra elements from the rotation matrix (3x3)
    """
    trace = numpy.matrix.trace(rotationMatrix)
    mod_w = math.acos((trace-1)/2)
    lie = mod_w * ((1/(2*math.sin(mod_w))) * numpy.matrix([[(rotationMatrix[2,1] - rotationMatrix[1,2])],[(rotationMatrix[0,2] - rotationMatrix[2,0])],[(rotationMatrix[1,0] - rotationMatrix[0,1])]]))
    return lie
    
def benjaminOlindeRodriguez(epsilon):
    """
    """
    u_hat = numpy.matrix([[0,-epsilon[2],epsilon[1]],[epsilon[2],0,-epsilon[0]],[-epsilon[1],epsilon[0],0]])
    u_mod = numpy.linalg.norm(epsilon)
    return numpy.transpose(numpy.eye(3) + ((numpy.transpose(u_hat)/u_mod)*math.sin(u_mod)) + (((numpy.transpose(u_hat)/u_mod)**2)*(1 - math.cos(u_mod))))
    

def checkRotationSanity(rotMatrix):
    ortho_identity = numpy.dot(rotMatrix,numpy.transpose(rotMatrix))
    ortho_identity = ortho_identity.round() + 0
    determinant = numpy.linalg.det(rotMatrix)

    if (((ortho_identity==numpy.eye(3)).all()) and (determinant.round()==1)):
        return 1
    else:
        return 0
    
    
    
    

rot = rotationMatrix(45,50,80,123)
print 'Initial rotation Matrix:'
print rot
lie = eulerLogarithm(rot)
print 'Lie: '
print lie
rot_back = benjaminOlindeRodriguez(lie)
print 'Re-traced rotation matrix:'
print rot_back
print 'Is rotation matrix valid? ',checkRotationSanity(rot_back)
