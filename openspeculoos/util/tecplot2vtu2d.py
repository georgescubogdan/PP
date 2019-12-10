#!/usr/bin/python

import sys

def locateZones(lines):
    zones = []
    for i in range(len(lines)):
        line = lines[i]
        if line[0:3]=="I =":
            lineNumber = i
            remainingLine = line[3:].partition(',')
            numX = int(remainingLine[0])
            remainingLine = remainingLine[2].partition('=')[2].partition(',')
            numY = int(remainingLine[0])
            zones.append( (lineNumber,numX,numY) )
    return zones

def writePoints(pointlist, fout):
    fout.write('<Points>\n')
    fout.write('<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n')
    for points in pointlist:
        fout.write(points[0] + ' ' + points[1] + ' ' + '0' + '\n')
    fout.write('</DataArray>\n');
    fout.write('</Points>\n');

def writeCells(numX, numY, fout):
    totElem = (numX-1)*(numY-1)
    fout.write('<Cells>\n')
    fout.write('<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n')
    for iY in range(numY-1):
        for iX in range(numX-1):
            p1 = iX + numX*iY
            p2 = p1+1
            p3 = p1 + numX + 1
            p4 = p3 - 1
            fout.write(str(p1) + ' ' + str(p2) + ' ' + str(p3) + ' ' + str(p4) + '\n')
    fout.write('</DataArray>\n');

    fout.write('<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n')
    for i in range(totElem):
        fout.write(str(4*(i+1))+' ')
    fout.write('\n</DataArray>\n');

    fout.write('<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n')
    for i in range(totElem):
        fout.write('9 ')
    fout.write('\n</DataArray>\n');
    fout.write('</Cells>\n')

def writeScalars(vx, vy, pressure, fout):
    fout.write('<PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n')

    fout.write('<DataArray type=\"Float64\" Name=\"Ux\" format=\"ascii\">\n')
    for ux in vx:
        fout.write(ux+' ')
    fout.write('\n</DataArray>\n')

    fout.write('<DataArray type=\"Float64\" Name=\"Uy\" format=\"ascii\">\n')
    for uy in vy:
        fout.write(uy+' ')
    fout.write('\n</DataArray>\n')

    fout.write('<DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n')
    for index in range(len(vx)):
        fout.write(vx[index]+' '+vy[index]+' '+'0'+'  ')
    fout.write('\n</DataArray>\n')

    fout.write('<DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">\n')
    for value in pressure:
        fout.write(value+' ')
    fout.write('\n</DataArray>\n')

    fout.write('</PointData>\n')


def formatZone(lines, numX, numY, fout):
    pointlist = []
    vXlist = []
    vYlist = []
    pressurelist = []
    pos = 0
    for y in range(numY):
        for x in range(numX):
            (pointX,pointY,vX,vY,pressure) = lines[pos].split()
            pointlist.append((pointX,pointY))
            vXlist.append(vX)
            vYlist.append(vY)
            pressurelist.append(pressure)
            pos = pos+1
    fout.write('<Piece NumberOfPoints=\"'+str(numX*numY)+'\" NumberOfCells=\"'+str((numX-1)*(numY-1))+'\">\n')
    writePoints(pointlist, fout)
    writeCells(numX, numY, fout)
    writeScalars(vXlist,vYlist,pressurelist, fout)
    fout.write('</Piece>\n\n')

if len(sys.argv)==2:
    fname_out = sys.argv[1] + ".vtu"
elif len(sys.argv)==3:
    fname_out = argv[2]
else:
    raise "Wrong number of arguments"


# MAIN
fin = open(sys.argv[1], 'r')
fout = open(fname_out, 'w')

lines = fin.read().splitlines()
zones = locateZones(lines)

fout.write('<?xml version=\"1.0\"?>\n')
fout.write('<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n')
fout.write('<UnstructuredGrid>\n\n')
for zone in zones:
    lineNumber = zone[0]
    numX = zone[1]
    numY = zone[2]
    formatZone(lines[lineNumber+1:lineNumber+1+numX*numY], numX, numY, fout)
fout.write('</UnstructuredGrid>\n')
fout.write('</VTKFile>\n')

