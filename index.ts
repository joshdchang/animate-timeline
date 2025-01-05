import { z } from "zod";
import { encode } from "@googlemaps/polyline-codec";
import sharp from "sharp";
import fs from "fs/promises";

// type defs
type Point = {
  timestamp: Date;
  coordinates: [number, number];
};
type Frame = {
  start: Date;
  end: Date;
  encodedLines: string[];
  center: [number, number];
  zoom: number;
};

// animation config
const animationStart = new Date("2024-01-01T00:00:00Z");
const animationEnd = new Date("2025-01-01T00:00:00Z");

const width = 1200;
const height = 800;

// get the data path from the command line args
const dataPath = Bun.argv[2];
if (!dataPath) {
  console.error("Usage: bun run plot.ts <location-history.json>");
  process.exit(1);
}

// get mapbox access token from env
const mapboxAccessToken = Bun.env.MAPBOX_ACCESS_TOKEN;
if (!mapboxAccessToken) {
  console.error("Please set the MAPBOX_ACCESS_TOKEN environment variable");
  process.exit(1);
}

// check that ffmpeg is installed
console.log("Checking for ffmpeg...");
const ffmpegInstalled = await Bun.$`ffmpeg -version`;
if (ffmpegInstalled.exitCode !== 0) {
  console.error("ffmpeg is not installed. Please install it and try again.");
  console.error("On macOS, you can install it with 'brew install ffmpeg'");
  process.exit(1);
}

// clear the frames directory
console.log("Clearing  data and preparing directories...");
await fs.rmdir("frames", { recursive: true });
await fs.mkdir("frames", { recursive: true });
await fs.rmdir("labelled_frames", { recursive: true });
await fs.mkdir("labelled_frames", { recursive: true });

// load the data
console.log("Loading location data from file...");
const rawData = await Bun.file(dataPath).json();
console.log(`Loaded ${rawData.length} location entries`);

// parse and validate the data sing Zod
console.log("Parsing and validating data...");
const locationSchema = z.array(
  z.object({
    endTime: z.string().transform((date) => new Date(date)),
    startTime: z.string().transform((date) => new Date(date)),
    timelinePath: z
      .array(
        z.object({
          point: z.string().transform((point) => {
            const [lat, lng] = point.split("geo:")[1].split(",");
            return [parseFloat(lat), parseFloat(lng)] as [number, number];
          }),
          durationMinutesOffsetFromStartTime: z.string().transform(Number),
        })
      )
      .optional(),
  })
);
const locationData = locationSchema.parse(rawData);

const points: Point[] = [];
locationData.forEach((location) => {
  location.timelinePath?.forEach((point) => {
    points.push({
      timestamp: new Date(
        location.startTime.getTime() +
          point.durationMinutesOffsetFromStartTime * 60 * 1000
      ),
      coordinates: point.point,
    });
  });
});
console.log(`Parsed and validated ${points.length} points`);

// interpolate the points
console.log(`Interpolating points...`);
const interpolatedPoints = generateInterpolatedPoints(points);
console.log(`Interpolated ${interpolatedPoints.length} points`);

// generate frame data
console.log("Generating frame data...");
const frames = generateFrames(interpolatedPoints, animationStart, animationEnd);
console.log(`Generated ${frames.length} frames`);

// smooth the zoom values
console.log("Smoothing zoom values...");
smoothZooms(frames, 168, 15.0);
console.log("Smoothed zoom values");

// build url for each frame
const accessToken =
  "pk.eyJ1Ijoiam9zaGNoYW5nMDQiLCJhIjoiY2p3c2c5NDdsMDEyOTQwcXhyc245eTYwcSJ9.CKHrw2ddbsyG9bdOlr0TvQ";
const stroke = "ff0000";
const weight = 3;

console.log("Generated image URLs");
const images = frames.map((frame, index) => {
  const [latCenter, lngCenter] = frame.center;
  return {
    url: `https://api.mapbox.com/styles/v1/mapbox/streets-v12/static/${frame.encodedLines
      .map(
        (encodedLine) => `path-${weight}+${stroke}(${encodeURIComponent(encodedLine)})`
      )
      .join(",")}/${lngCenter},${latCenter},${
      frame.zoom
    },0/${width}x${height}@2x?access_token=${accessToken}`,
    index,
  };
});
console.log("Generated image URLs");

// fetch & save
console.log("Fetching and saving frames...");
for (const { url, index } of images) {
  // console.log(`Fetching frame ${index} at ${url}`);
  fetch(url)
    .then((res) => res.arrayBuffer())
    .then((buffer) => {
      Bun.write(`./frames/frame-${index}.png`, buffer);
      console.log(`Frame ${index} saved`);
    })
    .catch((error) => {
      console.error(`Error fetching frame ${index}:`, error);
    });
}
console.log("All frames saved");

// add labels
console.log("Adding date labels to frames...");
const framesFiles = await fs.readdir("frames");

for (const frameFile of framesFiles) {
  const framePath = `frames/${frameFile}`;
  const frameImage = sharp(framePath);

  const frameNumber = parseInt(frameFile.split("_")[1].split(".")[0]);

  // every frame is one hour offset from the previous one, starting at animationStart
  const timestamp = animationStart.getTime() + frameNumber * 60 * 60 * 1000;
  const dateString = new Date(timestamp).toLocaleDateString("en-US", {
    timeZone: "UTC",
    month: "short",
    day: "numeric",
  });

  const frameWithText = frameImage.composite([
    {
      input: "box.png",
      gravity: "northwest",
    },
    {
      input: {
        text: {
          text: `<span background='white' foreground='black'>${dateString}</span>`,
          dpi: 300,
          rgba: true,
        },
      },
      blend: "darken",
      gravity: "northwest",
      left: 40,
      top: 30,
    },
  ]);

  frameWithText.toFile(`labelled_frames/${frameFile}`).then(() => {
    console.log(`Added label to frame ${frameNumber}`);
  });
}
console.log("All labels added");

// create the video
console.log("Creating video...");
Bun.$`ffmpeg -r 30 -f image2 -s ${width * 2}x${
  height * 2
} -start_number 1 -i labelled_frames/frame_%d.png -vframes ${
  frames.length - 1
} -vcodec libx264 -crf 25 -pix_fmt yuv420p ./movie.mp4`;
console.log("Video created");

// lib

function addMinutes(date: Date, minutes: number) {
  return new Date(date.getTime() + minutes * 60_000);
}

function interpolateGeodesic(
  lat1Deg: number,
  lng1Deg: number,
  lat2Deg: number,
  lng2Deg: number,
  t: number
): [number, number] {
  // convert degrees to radians
  const lat1 = (lat1Deg * Math.PI) / 180;
  const lng1 = (lng1Deg * Math.PI) / 180;
  const lat2 = (lat2Deg * Math.PI) / 180;
  const lng2 = (lng2Deg * Math.PI) / 180;

  // convert to cartesian coordinates on unit sphere
  const x1 = Math.cos(lat1) * Math.cos(lng1);
  const y1 = Math.cos(lat1) * Math.sin(lng1);
  const z1 = Math.sin(lat1);

  const x2 = Math.cos(lat2) * Math.cos(lng2);
  const y2 = Math.cos(lat2) * Math.sin(lng2);
  const z2 = Math.sin(lat2);

  // angle between the two vectors
  let dot = x1 * x2 + y1 * y2 + z1 * z2;
  // clamp
  if (dot < -1) dot = -1;
  if (dot > 1) dot = 1;

  const alpha = Math.acos(dot);

  // if alpha=0, the points are the same
  if (alpha === 0) {
    return [lat1Deg, lng1Deg];
  }

  const sinAlpha = Math.sin(alpha);
  const ratioA = Math.sin((1 - t) * alpha) / sinAlpha;
  const ratioB = Math.sin(t * alpha) / sinAlpha;

  // interpolated coordinates in 3D
  const x = ratioA * x1 + ratioB * x2;
  const y = ratioA * y1 + ratioB * y2;
  const z = ratioA * z1 + ratioB * z2;

  // convert back to spherical
  const interpLat = Math.atan2(z, Math.sqrt(x * x + y * y));
  const interpLng = Math.atan2(y, x);

  // convert radians back to degrees
  const latDeg = (interpLat * 180) / Math.PI;
  const lngDeg = (interpLng * 180) / Math.PI;

  return [latDeg, lngDeg];
}

function interpolateCoordinates(prev: Point, next: Point, time: Date): [number, number] {
  console.assert(
    prev.timestamp.getTime() <= next.timestamp.getTime(),
    `prev.timestamp (${prev.timestamp.toISOString()}) must be <= next.timestamp (${next.timestamp.toISOString()})`
  );

  const totalMs = next.timestamp.getTime() - prev.timestamp.getTime();
  const elapsedMs = time.getTime() - prev.timestamp.getTime();

  // edge case: if both timestamps are identical, just return prev's coords
  if (totalMs === 0) {
    return prev.coordinates;
  }

  console.assert(
    elapsedMs >= 0 && elapsedMs <= totalMs,
    `Time (${time.toISOString()}) is out of range for interpolation between ${prev.timestamp.toISOString()} and ${next.timestamp.toISOString()}`
  );

  const ratio = elapsedMs / totalMs;

  const [prevLat, prevLng] = prev.coordinates;
  const [nextLat, nextLng] = next.coordinates;

  return interpolateGeodesic(prevLat, prevLng, nextLat, nextLng, ratio);
}

function generateInterpolatedPoints(points: Point[]): Point[] {
  console.assert(points.length > 0, "Input array 'points' must not be empty");

  // sort the points by timestamp (if not already sorted)
  points.sort((a, b) => a.timestamp.getTime() - b.timestamp.getTime());

  // check that the points are sorted and have no duplicates
  for (let idx = 0; idx < points.length - 1; idx++) {
    console.assert(
      points[idx].timestamp.getTime() <= points[idx + 1].timestamp.getTime(),
      `Points are not sorted or have duplicates at index ${idx}`
    );
  }

  console.assert(points.length > 1, "Input array 'points' must have at least 2 points");

  const startTime = points[0].timestamp;
  const endTime = points[points.length - 1].timestamp;

  console.assert(
    startTime.getTime() <= endTime.getTime(),
    "Start time must be <= end time"
  );

  const interpolated: Point[] = [];

  // pointer for the prev point
  let i = 0;

  // walk minute-by-minute from start to end
  let currentTime = new Date(startTime);
  while (currentTime <= endTime) {
    // if the current time is beyond points[i+1], then move to the next pair
    // but make sure we don't run off the end of the array
    while (
      i < points.length - 1 &&
      currentTime.getTime() > points[i + 1].timestamp.getTime()
    ) {
      i++;
    }

    // if we're at the very last point, just use its coords
    if (i === points.length - 1) {
      interpolated.push({
        timestamp: new Date(currentTime),
        coordinates: points[i].coordinates,
      });
    } else {
      console.assert(i + 1 < points.length, "Index out of range for interpolation");

      // interpolate between points[i] (prev) and points[i+1] (next)
      const prevPoint = points[i];
      const nextPoint = points[i + 1];

      console.assert(
        currentTime.getTime() >= prevPoint.timestamp.getTime(),
        "currentTime cannot be before prevPoint.timestamp"
      );
      console.assert(
        currentTime.getTime() <= nextPoint.timestamp.getTime(),
        "currentTime cannot be after nextPoint.timestamp"
      );

      const coords = interpolateCoordinates(prevPoint, nextPoint, currentTime);

      interpolated.push({
        timestamp: new Date(currentTime),
        coordinates: coords,
      });
    }

    // advance time by 5 minutes
    currentTime = addMinutes(currentTime, 5);
  }

  return interpolated;
}

function weightedLatLngCenter(
  weightedPoints: { lat: number; lng: number; weight: number }[]
): [number, number] {
  let sumWeight = 0;
  let sumLatWeighted = 0;
  let sumLngWeighted = 0;

  for (const { lat, lng, weight } of weightedPoints) {
    sumLatWeighted += lat * weight;
    sumLngWeighted += lng * weight;
    sumWeight += weight;
  }

  if (sumWeight === 0) {
    return [0, 0]; // fallback
  }

  const centerLat = sumLatWeighted / sumWeight;
  const centerLng = sumLngWeighted / sumWeight;
  return [centerLat, centerLng];
}

function computeWeightedLatLngStdDevKm(
  weightedPoints: { lat: number; lng: number; weight: number }[],
  center: [number, number]
): { latStdKm: number; lngStdKm: number } {
  const [centerLatDeg, centerLngDeg] = center;

  const DEG_TO_KM_LAT = 111; // approx 111 km per degree of latitude
  const degToKmLng = 111 * Math.cos((centerLatDeg * Math.PI) / 180);

  let sumWeight = 0;
  let sumLatSqWeighted = 0;
  let sumLngSqWeighted = 0;

  for (const { lat, lng, weight } of weightedPoints) {
    const dLatDeg = lat - centerLatDeg;
    const dLngDeg = lng - centerLngDeg;

    const dLatKm = dLatDeg * DEG_TO_KM_LAT;
    const dLngKm = dLngDeg * degToKmLng;

    sumLatSqWeighted += weight * (dLatKm * dLatKm);
    sumLngSqWeighted += weight * (dLngKm * dLngKm);
    sumWeight += weight;
  }

  if (sumWeight === 0) {
    return { latStdKm: 0, lngStdKm: 0 };
  }

  const latVariance = sumLatSqWeighted / sumWeight;
  const lngVariance = sumLngSqWeighted / sumWeight;

  // std dev
  const latStdKm = Math.sqrt(latVariance);
  const lngStdKm = Math.sqrt(lngVariance);

  return { latStdKm, lngStdKm };
}

function computeZoomFromLatLngStdDevs(
  weightedPoints: { lat: number; lng: number; weight: number }[]
): number {
  const center = weightedLatLngCenter(weightedPoints);
  const { latStdKm, lngStdKm } = computeWeightedLatLngStdDevKm(weightedPoints, center);
  const aspectRatio = width / height;
  const maxStdKm = Math.max(latStdKm * aspectRatio, lngStdKm);
  return calculateZoomFromDistance(maxStdKm);
}

function calculateZoomFromDistance(distKm: number): number {
  if (distKm <= 0) {
    return 6; // default
  }

  // some random math and magic
  const desiredSpanKm = 4 * distKm;
  const fraction = 40000 / desiredSpanKm;
  const pxRatio = width / 256;
  let zoom = Math.log2(fraction * pxRatio);

  // clamp & tweak
  zoom = Math.max(3, Math.min(11.8, zoom / 1.4 + 0.4));
  return zoom;
}

function weightedGeodesicCenter(
  weightedPoints: { lat: number; lng: number; weight: number }[]
): [number, number] {
  let xSum = 0;
  let ySum = 0;
  let zSum = 0;
  let totalWeight = 0;

  for (const { lat, lng, weight } of weightedPoints) {
    // to radians
    const latRad = (lat * Math.PI) / 180;
    const lngRad = (lng * Math.PI) / 180;

    // Cartesian on unit sphere
    const cosLat = Math.cos(latRad);
    const x = cosLat * Math.cos(lngRad);
    const y = cosLat * Math.sin(lngRad);
    const z = Math.sin(latRad);

    // accumulate weighted sums
    xSum += x * weight;
    ySum += y * weight;
    zSum += z * weight;
    totalWeight += weight;
  }

  // edge case: all zero weights
  if (totalWeight === 0) {
    return [0, 0];
  }

  // avg
  xSum /= totalWeight;
  ySum /= totalWeight;
  zSum /= totalWeight;

  // back to lat/lng
  const hyp = Math.sqrt(xSum * xSum + ySum * ySum);
  const centerLat = Math.atan2(zSum, hyp) * (180 / Math.PI);
  const centerLng = Math.atan2(ySum, xSum) * (180 / Math.PI);

  return [centerLat, centerLng];
}

function breakAntimeridianSimple(coords: [number, number][]): Array<[number, number][]> {
  if (coords.length <= 1) {
    return [coords];
  }

  const segments: Array<[number, number][]> = [];
  let currentSegment: [number, number][] = [coords[0]];

  for (let i = 1; i < coords.length; i++) {
    const [latPrev, lngPrev] = coords[i - 1];
    const [latCurr, lngCurr] = coords[i];

    // if we jump from >175 to <−175, or <−175 to >175, consider it a crossing, nice and simple
    if ((lngPrev > 175 && lngCurr < -175) || (lngPrev < -175 && lngCurr > 175)) {
      // end the current segment
      segments.push(currentSegment);

      // start a new segment
      currentSegment = [[latCurr, lngCurr]];
    } else {
      // add the next point to the current segment
      currentSegment.push([latCurr, lngCurr]);
    }
  }

  // push the last segment
  if (currentSegment.length > 0) {
    segments.push(currentSegment);
  }

  return segments;
}

function gaussianWeight(timestampMs: number, centerMs: number, halfWindow: number) {
  const dist = Math.abs(timestampMs - centerMs);
  // Pick some sigma ~ (halfWindow / 2) or whichever you prefer
  const sigma = halfWindow / 2;
  return Math.exp(-0.5 * (dist / sigma) ** 2);
}

function generateFrames(
  points: Point[],
  animationStart: Date,
  animationEnd: Date,
  stepInHours: number = 1,
  pathWindowInHours: number = 24,
  centerWindowInHours: number = 120
) {
  // prefix sums for faster path-encoding
  const latPrefix = new Array<number>(points.length + 1).fill(0);
  const lngPrefix = new Array<number>(points.length + 1).fill(0);
  for (let i = 0; i < points.length; i++) {
    latPrefix[i + 1] = latPrefix[i] + points[i].coordinates[0];
    lngPrefix[i + 1] = lngPrefix[i] + points[i].coordinates[1];
  }

  const frames: Array<Frame> = [];

  // for slicing the path (24h window by default)
  let pathLeft = 0;
  let pathRight = 0;

  // for slicing the center (120h window by default)
  let centerLeft = 0;
  let centerRight = 0;

  const pathWindowSide = (pathWindowInHours / 2) * 60 * 60_000;
  const centerWindowSide = (centerWindowInHours / 2) * 60 * 60_000;

  let i = 0;
  while (true) {
    const frameMidMs = animationStart.getTime() + i * stepInHours * 60 * 60_000;
    if (frameMidMs > animationEnd.getTime()) {
      break;
    }

    const pathStartMs = frameMidMs - pathWindowSide;
    const pathEndMs = frameMidMs + pathWindowSide;

    // increment until where we want to be
    while (
      pathLeft < points.length &&
      points[pathLeft].timestamp.getTime() < pathStartMs
    ) {
      pathLeft++;
    }
    while (
      pathRight < points.length &&
      points[pathRight].timestamp.getTime() <= pathEndMs
    ) {
      pathRight++;
    }

    const pathPoints = points.slice(pathLeft, pathRight);
    const pathCoords = pathPoints.map((p) => p.coordinates);

    const segments = breakAntimeridianSimple(pathCoords);

    // encode each segment separately
    const encodedLines = segments.map((segment) => encode(segment));

    // weighted center calculation
    const centerStartMs = frameMidMs - centerWindowSide;
    const centerEndMs = frameMidMs + centerWindowSide;

    // increment again until where we want to be
    while (
      centerLeft < points.length &&
      points[centerLeft].timestamp.getTime() < centerStartMs
    ) {
      centerLeft++;
    }
    while (
      centerRight < points.length &&
      points[centerRight].timestamp.getTime() <= centerEndMs
    ) {
      centerRight++;
    }

    // weighted sums in 3d
    const weightedPoints = [];

    for (let idx = centerLeft; idx < centerRight; idx++) {
      const { timestamp, coordinates } = points[idx];
      const [lat, lng] = coordinates;
      const w = gaussianWeight(timestamp.getTime(), frameMidMs, centerWindowSide);
      if (w > 0) {
        weightedPoints.push({ lat, lng, weight: w });
      }
    }

    const center: [number, number] = weightedGeodesicCenter(weightedPoints);

    const zoom = computeZoomFromLatLngStdDevs(weightedPoints);

    frames.push({
      start: new Date(pathStartMs),
      end: new Date(pathEndMs),
      encodedLines,
      center,
      zoom,
    });

    i++;
  }

  return frames;
}

function smoothZooms(frames: Frame[], radius = 2, sigma = 1.0) {
  // generate the kernel weights
  const kernelSize = radius * 2 + 1;
  const kernel = new Array<number>(kernelSize);

  let sum = 0;
  for (let i = -radius; i <= radius; i++) {
    const x = (i * i) / (2 * sigma * sigma);
    // basic Gaussian kernel: e^{-x}
    const w = Math.exp(-x);
    kernel[i + radius] = w;
    sum += w;
  }

  // normalize the kernel so that all weights sum to 1
  for (let i = 0; i < kernelSize; i++) {
    kernel[i] /= sum;
  }

  // for each frame, apply the kernel
  const smoothedZooms = new Array<number>(frames.length);
  for (let i = 0; i < frames.length; i++) {
    let weightedSum = 0;
    for (let j = -radius; j <= radius; j++) {
      const neighborIndex = i + j;
      if (neighborIndex >= 0 && neighborIndex < frames.length) {
        // kernelIndex is the offset within the kernel array
        const kernelIndex = j + radius;
        weightedSum += frames[neighborIndex].zoom * kernel[kernelIndex];
      }
    }

    smoothedZooms[i] = weightedSum;
  }

  // store the smoothed zoom values back into frames
  for (let i = 0; i < frames.length; i++) {
    frames[i].zoom = smoothedZooms[i];
  }
}
