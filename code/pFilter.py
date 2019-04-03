import numpy as np
import math
import sys
import cv2
import pdb

# /**
#  * Main function of the propagation filter
#  * @param RPtr pointer to an array containing reference matrix 
#  * @param APtr pointer to an array containing input matrix 
#  * @param w window radius
#  * @param sigma_d
#  * @param sigma_r
#  * @param nRows the number of rows of R and A
#  * @param nCols the number of columns of R and A
#  * @param zR the number of channels of R
#  * @param zA the number of channels of A
#  */

def pfilter_matrix(RPtr, APtr, w, sigma_d, sigma_r, nRows,  nCols, zR, zA) :
    
    outImg = np.zeros((nRows, nCols, zA)) 
    for i in range(nRows):      
        for j in range(nCols):
            point = [i, j]             
            resultPixelPtr = calculateWeightedAverage(RPtr, APtr, w, sigma_d, sigma_r, point, nRows, nCols, zR, zA);
            outImg[i][j] = resultPixelPtr
            # pdb.set_trace()

    return outImg


# /**
#  * Calculate weighted average of pixels surrounding the point
#  * @param resultPixelPtr the array (with the size = zA) contains the result pixel values 
#  * @param RPtr
#  * @param APtr
#  * @param w
#  * @param sigma_d
#  * @param sigma_r
#  * @param point
#  * @param m the number of rows of R and A
#  * @param n the number of columns of R and A
#  * @return void
#  */
def calculateWeightedAverage(RPtr, APtr, w, sigma_d, sigma_r, point, m, n, zR, zA) :

    resultPixelPtr = np.zeros(zA)
    wSize = 2 * w + 1;
    logWeightWindowPtr = np.zeros((wSize, wSize))
    totalWeight = 0
    totalWeightedSum = np.zeros(zA)
    
    # Calculate the weight of the Center Point
    logWeightWindowPtr[w][w] = 0
    totalWeight = 1.0
    totalWeightedSum = APtr[point[0]][point[1]] 

    # // Calculate from distance 1 to w
    for r in range(1, w+1):
        for dp in range(r+1):
            for pSign in range(-1,2,2): # sign = -1, 1
                p = pSign*dp;
                for qSign in range(-1,2,2): # sign = -1, 1
                    q = qSign * (r - dp)
                    # // check boundary
                    if (point[0] + p < 0 or point[0] + p > m - 1 or point[1] + q < 0 or point[1] + q > n - 1) :
                        continue

                    # decide fromLocal (the parent pixel t-1)
                    fromLocal = [0, 0]
                    if (p * q == 0): # on the x or y axis
                        if (p == 0) :
                            fromLocal[0] = p
                            fromLocal[1] = q - qSign
                        else : # q == 0
                            fromLocal[0] = p - pSign;
                            fromLocal[1] = q;

                    else : # p*q != 0 (other pixels)
                        # if r is odd -> p , else -> q
                        if (r % 2 != 0) :
                            fromLocal[0] = p;
                            fromLocal[1] = q - qSign;
                        else :
                            fromLocal[0] = p - pSign;
                            fromLocal[1] = q;

                    # calculate log weight
                    toLocal = np.array([p, q]);
                    logWeight = calculateLogWeight(RPtr, sigma_d, sigma_r, logWeightWindowPtr, w, point, fromLocal, toLocal);
                    logWeightWindowPtr[w+p][w+q] = logWeight
                    weight = math.exp(logWeight);               
                    
                    totalWeight = totalWeight + weight;
                    totalWeightedSum = totalWeightedSum + weight*APtr[point[0]+p][point[1]+q] 
                    
                    # ensure pixels on the axis is calculated only one time 
                    if q == 0 :
                        break
                    
                # ensure pixels on the axis is calculated only one time 
                if p == 0 :
                    break
    
    
    # Calculate result pixel value
    resultPixelPtr = totalWeightedSum / totalWeight; 

    return resultPixelPtr


# /**
#  * Calculate the log relationship of two points in R
#  * @param RPtr
#  * @param sigma
#  * @param fromPoint
#  * @param toPoint
#  * @return 
#  */
def calculateLogRelationship(RPtr, sigma, fromPoint, toPoint)  :

    # pdb.set_trace()
    diff = RPtr[fromPoint[0]][fromPoint[1]] - RPtr[toPoint[0]][toPoint[1]]
    distanceSquare = np.sum(np.power(diff, 2))

    return -1 * distanceSquare / (2 * math.pow(sigma, 2))


# /**
#  * Calculate the log-weight between centerPoint and the toPoint through fromPoint
#  * @param RPtr
#  * @param sigma_d
#  * @param sigma_r
#  * @param fromPoint
#  * @param toPoint
#  * @return 
#  */
def calculateLogWeight(RPtr, sigma_d, sigma_r, logWeightWindowPtr, w, centerPoint, fromLocal, toLocal) :

    fromPoint = [centerPoint[0] + fromLocal[0], centerPoint[1] + fromLocal[1]]
    toPoint = [centerPoint[0] + toLocal[0], centerPoint[1] + toLocal[1]]

    pathLogProb = calculateLogRelationship(RPtr, sigma_d, fromPoint, toPoint)
    rangeLogProb = calculateLogRelationship(RPtr, sigma_r, centerPoint, toPoint)
    logWeight = logWeightWindowPtr[fromLocal[0] + w][fromLocal[1] + w] + pathLogProb + rangeLogProb

    return logWeight


def main(argv):
    originalImage = cv2.imread(argv[1])
    referenceImage = cv2.imread(argv[2])
    window_radius = int(argv[3])
    sigma_d = int(argv[4])
    sigma_r = int(argv[5])
    nRows, nCols, zA = originalImage.shape
    print nRows, nCols, zA
    zR = referenceImage.shape[-1]
    finalImage = pfilter_matrix(referenceImage, originalImage, window_radius, sigma_d, sigma_r, nRows,  nCols, zR, zA)
    return finalImage


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print "Usage : python pFilter.py <Original Image> <Reference Image> <Window Radius> <Sigma_d> <Sigma_r>"
    else :
        imgName = sys.argv[1].split('/')[-1]
        imgName = imgName.split('.')[0] + '_' + sys.argv[3] + '_' + sys.argv[4] + '_'
        print imgName
        finalImage = main(sys.argv)
        cv2.imwrite(imgName+'01.png',finalImage)
        # pdb.set_trace()
        finalImage = finalImage.astype(int)
        cv2.imwrite(imgName + '02.png',finalImage)
        # cv2.imshow('image',finalImage)
        # cv2.waitKey(0)
        # cv2.destroyAllWindows()
