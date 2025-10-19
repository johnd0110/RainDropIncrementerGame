import random
import math
import time
from argparse import ArgumentError
from tkinter import *
from tkinter import ttk

def convertListOfTuplesToList(listOfTuples):
    return [item for tup in listOfTuples for item in tup]

def getNumSign(num):
    return 0 if num == 0 else math.copysign(1, num)

class rectangleWalker:
    minimumWalkDistance = 20
    def __init__(self, canvasWidth, canvasHeight):
        self.canvasWidth = canvasWidth
        self.canvasHeight = canvasHeight
        self.quarterCanvasWidth = self.canvasWidth // 4
        self.quarterCanvasHeight = self.canvasHeight // 4

        # TODO: Fix the rectangle so it isn't as tall
        # Top Left Vertice will start out in a random spot in the right half of the top left quadrant
        startingRectangleCoord = (random.randint(0, self.canvasWidth), random.randint(self.quarterCanvasHeight + self.quarterCanvasHeight//2, self.canvasHeight//2))
        # Right side X coordinate will be in the left half of the top right quadrant
        nextRectangleXCoord = startingRectangleCoord[0] + random.randint(self.quarterCanvasWidth//2, self.quarterCanvasWidth)
        # Bottom Y coordinate will be in the upper half of the bottom 2 quadrants
        nextRectangleYCoord = startingRectangleCoord[1] + random.randint(self.quarterCanvasHeight//4, self.quarterCanvasHeight//2)

        # Coordinates in order of top left -> top right -> Bottom Right -> Bottom Left
        self.rectangleCoordinates = [startingRectangleCoord,
                                    (nextRectangleXCoord, startingRectangleCoord[1]),
                                    (nextRectangleXCoord, nextRectangleYCoord),
                                    (startingRectangleCoord[0], nextRectangleYCoord),
        ]
        self.index = 0
        self.currentCoordinate = self.rectangleCoordinates[self.index]

    def __iter__(self):
        return self

    def __next__(self):
        # TODO: Fix the walk so that it doesn't have the chance to walk the entire length of the rectangle side in one go
        def progressWalk(directionX, directionY):
            if (directionX or directionY) not in (-1, 0, 1):
                raise ArgumentError("Direction X/Y must be a -1, 0, or 1.")
            return (self.currentCoordinate[0] + (0 if directionX == 0 else random.randint(rectangleWalker.minimumWalkDistance, self.quarterCanvasWidth//4) * directionX),
                    self.currentCoordinate[1] + (0 if directionY == 0 else random.randint(rectangleWalker.minimumWalkDistance, self.canvasHeight//4) * directionY))

        # Movement is strictly increasing (+X or +Y)
        nextIndex = self.index + 1

        if nextIndex < 4:
            nextVertice = self.rectangleCoordinates[nextIndex]
        elif nextIndex == 4:
            nextVertice = self.rectangleCoordinates[0]
        else:
            raise StopIteration

        # TODO: Fix walk to avoid issue of a given walk being too short (like 1 pixel).
        #  Can occur when the prior coordinate was very close to the next vertice and the natural progression caused the new coordinate to be clamped to the next vertice
        # We have moved past the next vertice, need to clamp back to the next vertice and increment the index
        if (nextIndex < 3 and self.currentCoordinate > nextVertice) or (nextIndex >= 3 and self.currentCoordinate < nextVertice):
            self.currentCoordinate = nextVertice
            self.index += 1
        # The current coordinate is the next vertice, thus we can just increment the index and return the next vertice as the current coordinate
        elif self.currentCoordinate == nextVertice:
            self.index += 1

        coordinate = self.currentCoordinate
        # The Index might have been incremented based on earlier logic
        currentNextIndex = self.index + 1
        match currentNextIndex:
            case 1:
                # Movement is in the positive X Direction
                self.currentCoordinate = progressWalk(1, 0)
            case 2:
                # Movement is in the positive Y direction
                self.currentCoordinate = progressWalk(0, 1)
            case 3:
                # Movement is in the negative X direction
                self.currentCoordinate = progressWalk(-1, 0)
            case 4:
                # Movement is in the negative Y direction
                self.currentCoordinate = progressWalk(0, -1)
            case 5:
                pass
            case _:
                raise ValueError(f"Index: {currentNextIndex} should be less than 5")

        return coordinate

class rainCloud:
    cloudTagPrefix = "cloud"
    maxCloudTagIds = 100
    thunderTagPrefix = "thunder"
    maxThunderTagIds = 100
    driftMovementAmount = 50
    driftMovementMaxCount = 15
    driftMovementDelay = 500 #In ms

    def __init__(self, canvas):
        self.dirty = False
        self.canvas = canvas
        self.cloudTagId = next(rainCloudTagGenerator)
        self.rainDropDict = {}

        canvasSizeWidth, canvasSizeHeight = self.retrieveCanvasSize()
        rectWalk = rectangleWalker(canvasSizeWidth, canvasSizeHeight)
        arcCoords = [coord for coord in rectWalk]

        self.cloudRectId = self.canvas.create_rectangle(convertListOfTuplesToList([rectWalk.rectangleCoordinates[0], rectWalk.rectangleCoordinates[2]]), fill="gray", tags=(self.cloudTagId))

        for ind, coordinate in enumerate(arcCoords):
            nextIndex = ind + 1
            if nextIndex >= len(arcCoords):
                break

            arcCoordPair = [coordinate, arcCoords[nextIndex]]

            walkDirection = tuple(getNumSign(zippedPair[1] - zippedPair[0]) for zippedPair in zip(arcCoordPair[0], arcCoordPair[1]))

            arcHeight = random.randint(rainDropIncrementer.minimumArcHeight, rainDropIncrementer.maximumArcHeight) // 2
            match walkDirection:
                case (1, 0):
                    arcCoordPair[0] = (arcCoordPair[0][0], arcCoordPair[0][1] - arcHeight)
                    arcCoordPair[1] = (arcCoordPair[1][0], arcCoordPair[1][1] + arcHeight)
                    arcCoordPairList = convertListOfTuplesToList(arcCoordPair)
                    self.canvas.create_oval(arcCoordPairList, fill="gray", outline="", tags=(self.cloudTagId))
                    self.canvas.create_arc(arcCoordPairList, style=ARC, outline="black", start=0, extent=180, tags=(self.cloudTagId))
                case (0, 1):
                    arcCoordPair[0] = (arcCoordPair[0][0] - arcHeight, arcCoordPair[0][1])
                    arcCoordPair[1] = (arcCoordPair[1][0] + arcHeight, arcCoordPair[1][1])
                    arcCoordPairList = convertListOfTuplesToList(arcCoordPair)
                    self.canvas.create_oval(arcCoordPairList, fill="gray", outline="", tags=(self.cloudTagId))
                    self.canvas.create_arc(arcCoordPairList, style=ARC, outline="black", start=270, extent=180, tags=(self.cloudTagId))
                case (-1, 0):
                    arcCoordPair[0] = (arcCoordPair[0][0], arcCoordPair[0][1] - arcHeight)
                    arcCoordPair[1] = (arcCoordPair[1][0], arcCoordPair[1][1] + arcHeight)
                    arcCoordPairList = convertListOfTuplesToList(arcCoordPair)
                    self.canvas.create_oval(arcCoordPairList, fill="gray", outline="", tags=(self.cloudTagId))
                    self.canvas.create_arc(arcCoordPairList, style=ARC, outline="black", start=180, extent=180, tags=(self.cloudTagId))
                case (0, -1):
                    arcCoordPair[0] = (arcCoordPair[0][0] - arcHeight, arcCoordPair[0][1])
                    arcCoordPair[1] = (arcCoordPair[1][0] + arcHeight, arcCoordPair[1][1])
                    arcCoordPairList = convertListOfTuplesToList(arcCoordPair)
                    self.canvas.create_oval(arcCoordPairList, fill="gray", outline="", tags=(self.cloudTagId))
                    self.canvas.create_arc(arcCoordPairList, style=ARC, outline="black", start=90, extent=180, tags=(self.cloudTagId))
                case _:
                    raise ValueError(f"{walkDirection} is not a valid walking direction")

        self.animateDrift()

    def retrieveCanvasSize(self):
        actualCanvasSize = (self.canvas.winfo_width(), self.canvas.winfo_height())
        return actualCanvasSize

    @staticmethod
    def tagGenerator(tagPrefix, maxTagCount):
        existingTags = []

        def newTag(generateNewTag):
            if not generateNewTag: return None
            tagIds = [int(tagId.removeprefix(tagPrefix)) for tagId in existingTags]
            tag = tagPrefix + str(random.choice([i for i in range(0, maxTagCount) if i not in tagIds]))
            existingTags.append(tag)
            return tag

        tagId = None
        while True:
            generateNewTag = True
            if tagId in existingTags:
                #print("cleaned!")
                existingTags.remove(tagId)
                generateNewTag = False
            #print(f"Before Generating new cloud tag. Existing Cloud Tags: {existingCloudTags}")
            tagId = yield newTag(generateNewTag)
            #print(f"After Generating new cloud tag. Existing Cloud Tags: {existingCloudTags}")

    def animateDrift(self):
        cloudBoundingBox = self.canvas.bbox(self.cloudTagId)
        if cloudBoundingBox == "":
            raise ValueError("No Cloud Bounding Box found.")

        # When the cloud has reached the horizontal ends of the canvas, we mark them to be cleaned up and terminate the drifting animation
        if cloudBoundingBox[0] > self.canvas.winfo_width() or cloudBoundingBox[0] < 0:
            self.dirty = True
            return
        else:
            self.canvas.move(self.cloudTagId, random.randint(0, rainCloud.driftMovementAmount) * random.randrange(-1, 2, 1), 0)

        self.canvas.after(random.randint(rainCloud.driftMovementDelay//2, rainCloud.driftMovementDelay), self.animateDrift)

    def retrieveRandomCoordinateAlongCloudBottom(self):
        cloudRectCoords = self.canvas.coords(self.cloudRectId)
        return (random.randint(int(cloudRectCoords[0]), int(cloudRectCoords[2])), int(cloudRectCoords[-1]))

    def rainDropAnimation(self, rainDropId):
        if rainDropId is None:
            # If the cloud has been set to be cleaned up, don't create any more rain drops.
            if self.dirty:
                return
            # Create a rain drop and spawn it inside our rain cloud
            try:
                lineXCoord, lineYCoord = self.retrieveRandomCoordinateAlongCloudBottom()
            except IndexError:
                return
            newRainDropId = self.canvas.create_line(lineXCoord, lineYCoord, lineXCoord, lineYCoord + rainDropIncrementer.rainDropHeight, fill="blue")

            if newRainDropId in self.rainDropDict:
                raise KeyError(f"Rain Drop ID: {newRainDropId} somehow already exists for cloud tag ID: {self.cloudTagId}")
            # Created a new rain drop for our rain cloud so we need to initialize the movement count
            self.rainDropDict[newRainDropId] = 0
            self.canvas.after(rainDropIncrementer.rainDropAnimationDelay, self.rainDropAnimation, newRainDropId)
        else:
            # Processing next animation step of an existing rain drop
            if self.rainDropDict[rainDropId] < rainDropIncrementer.rainDropMaxMovementCount:
                # Move rain drop a step lower to simulate the falling animation
                self.canvas.move(rainDropId, 0, rainDropIncrementer.rainDropHeight)

                # To save some processing, we terminate the animation if the rain drop has fallen past the bottom of the canvas
                if self.canvas.coords(rainDropId)[-1] > self.canvas.winfo_height():
                    self.rainDropDict[rainDropId] = rainDropIncrementer.rainDropMaxMovementCount
                else:
                    self.rainDropDict[rainDropId] += 1

                # Schedule up the next animation step
                self.canvas.after(rainDropIncrementer.rainDropAnimationDelay, self.rainDropAnimation, rainDropId)
            else:
                # Animation is over, so we clean up the rain drop
                self.canvas.delete(rainDropId)
                _ = self.rainDropDict.pop(rainDropId)

    def drawThunder(self):
        if self.dirty: return
        startingXCoord, startingYCoord = self.retrieveRandomCoordinateAlongCloudBottom()

        def drawUpsideDownTrapezoid(bottomWidth, topWidth, height, color, tag):
            if topWidth >= bottomWidth:
                raise ArgumentError("Top Side Width >= bottom side width, cannot draw trapezoid.")
            xAdjustment = (bottomWidth - topWidth) // 2
            return self.canvas.create_polygon(0,0, bottomWidth,0, bottomWidth-xAdjustment, height, xAdjustment, height, fill=color, outline="", tags=(tag))

        def drawUpsideDownTriangle(baseWidth, height, color, tag):
            return self.canvas.create_polygon(0,0, baseWidth,0, baseWidth//2, height, fill=color, outline="", tags=(tag))

        thunderTag = next(thunderTagGenerator)

        trapezoidBottomWidth = 30
        trapezoidTopWidth = 20
        trapezoidHeight = 30

        _ = drawUpsideDownTrapezoid(trapezoidBottomWidth, trapezoidTopWidth, trapezoidHeight, "goldenrod", thunderTag)
        triangleID = drawUpsideDownTriangle(trapezoidTopWidth, trapezoidHeight*2, "goldenrod", thunderTag)
        self.canvas.move(triangleID, (trapezoidBottomWidth - trapezoidTopWidth), trapezoidHeight)
        self.canvas.moveto(thunderTag, startingXCoord, startingYCoord)

        self.canvas.after(500, self.cleanupThunder, thunderTag)

    def cleanupThunder(self, thunderTag):
        self.canvas.delete(thunderTag)
        thunderTagGenerator.send(thunderTag)

rainCloudTagGenerator = rainCloud.tagGenerator(rainCloud.cloudTagPrefix, rainCloud.maxCloudTagIds)
thunderTagGenerator = rainCloud.tagGenerator(rainCloud.cloudTagPrefix, rainCloud.maxCloudTagIds)

class rainDropIncrementer:
    startingColumn = 0
    startingRow = 0
    minimumArcHeight = 25
    maximumArcHeight = 50
    rainDropHeight = 30
    rainDropMaxMovementCount = 15
    rainDropAnimationDelay = 125 # In ms
    maxNumOfClouds = 4
    rainDropThunderThreshold = 5
    def __init__(self, tkRoot):
        self.root = tkRoot
        self.root.title("Rain Drop")

        self.mainframe = ttk.Frame(self.root, padding=(3, 3, 12, 12))
        self.mainframe.grid(column=rainDropIncrementer.startingColumn, row=rainDropIncrementer.startingRow, sticky=(N, W, E, S))

        self.cloudCanvas = Canvas(self.mainframe, width=500, height=500)
        self.cloudCanvas.grid(column=rainDropIncrementer.startingColumn+1,row=rainDropIncrementer.startingRow+1, sticky=(W,E), columnspan=4)

        self.cloudCanvas.update()
        self.clouds = {}
        self.cleanupAndGenerateClouds()

        ttk.Label(self.mainframe, text="Rain Drops:").grid(column=rainDropIncrementer.startingColumn+1, row=rainDropIncrementer.startingRow+2, sticky=(E))
        self.rainDropCount = IntVar()
        self.rainDropCount.set(0)
        ttk.Label(self.mainframe, textvariable=self.rainDropCount).grid(column=rainDropIncrementer.startingColumn+2, row=rainDropIncrementer.startingRow+2, sticky=(W))

        self.rainDropRate = DoubleVar()
        self.rainDropRate.set(0)
        ttk.Label(self.mainframe, textvariable=self.rainDropRate).grid(column=rainDropIncrementer.startingColumn+3, row=rainDropIncrementer.startingRow+2, sticky=(E))
        ttk.Label(self.mainframe, text="Rain Drops/s").grid(column=rainDropIncrementer.startingColumn+4, row=rainDropIncrementer.startingRow+2, sticky=(W))
        self.CalculateAndSetRainDropRate()
        self.generateThunder()

        ttk.Button(self.mainframe, text="Increment", command=self.Incrementer).grid(column=rainDropIncrementer.startingColumn+1, row=rainDropIncrementer.startingRow+3, sticky=(E, W), columnspan=4)

        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.mainframe.columnconfigure(1, weight=1)
        self.mainframe.columnconfigure(2, weight=1)
        for child in self.mainframe.winfo_children():
            child.grid_configure(padx=5, pady=5)

        self.root.bind('<Return>', self.Incrementer)

    def cleanupAndGenerateClouds(self):
        cloudsToCleanup = []
        for cloudTagId in self.clouds:
            cloud = self.getCloudForCloudTagId(cloudTagId)
            if cloud.dirty and len(cloud.rainDropDict) == 0:
                cloudsToCleanup.append(cloudTagId)

        for cloudTagId in cloudsToCleanup:
            self.cloudCanvas.delete(cloudTagId)
            self.clouds.pop(cloudTagId)
            rainCloudTagGenerator.send(cloudTagId)

        if len(self.clouds) < rainDropIncrementer.maxNumOfClouds:
            self.generateNewClouds()

        self.cloudCanvas.after(250, self.cleanupAndGenerateClouds)

    def generateNewClouds(self):
        maxAmountofCloudsToGenerate = rainDropIncrementer.maxNumOfClouds - len(self.clouds)

        if maxAmountofCloudsToGenerate < 0:
            raise ValueError("Negative amount of clouds to generate.")

        for i in range(0, random.randint(0, maxAmountofCloudsToGenerate)):
            self.addCloud()

    def addCloud(self):
        aRainCloud = rainCloud(self.cloudCanvas)
        self.clouds[aRainCloud.cloudTagId] = aRainCloud

    def getRandomCleanCloud(self):
        return random.choice([cloud for cloud in self.clouds.values() if not cloud.dirty])

    def getCloudForCloudTagId(self, cloudTagId):
        return self.clouds[cloudTagId]

    def CalculateAndSetRainDropRate(self, previousRainDropAmount=0, previousTimestamp=None):
        currentTimeStamp = time.time()
        currentRainDropAmount = self.rainDropCount.get()
        if previousTimestamp is not None:
            self.rainDropRate.set(round((currentRainDropAmount - previousRainDropAmount) / (currentTimeStamp - previousTimestamp), 2))

        self.root.after(1000, self.CalculateAndSetRainDropRate, currentRainDropAmount, currentTimeStamp)

    def generateThunder(self):
        if self.rainDropRate.get() >= rainDropIncrementer.rainDropThunderThreshold:
            self.getRandomCleanCloud().drawThunder()
        self.root.after(1000, self.generateThunder)

    def Incrementer(self, *args):
        try:
            self.rainDropCount.set(self.rainDropCount.get() + 1)
            self.getRandomCleanCloud().rainDropAnimation(None)
        except ValueError:
            pass

if __name__ == '__main__':
    root = Tk()
    rainDropIncrementer(root)
    root.mainloop()