import { z } from "zod";
import { encode } from "@googlemaps/polyline-codec";
import fs from "fs/promises";
import { animationEnd, animationStart } from "./config";


const dataPath = Bun.argv[2];
if (!dataPath) {
  console.error("Usage: bun run plot.ts <data.json>");
  process.exit(1);
}

await fs.mkdir("frames", { recursive: true });

const rawData = await Bun.file(dataPath).json();

const width = 1200;
const height = 800;

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

type Point = {
  timestamp: Date;
  coordinates: [number, number];
};

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

// get points for each minute from start to end by interpolating between the points

/**
 * Helper to add minutes to a Date.
 */
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
  // Convert degrees to radians
  const lat1 = (lat1Deg * Math.PI) / 180;
  const lng1 = (lng1Deg * Math.PI) / 180;
  const lat2 = (lat2Deg * Math.PI) / 180;
  const lng2 = (lng2Deg * Math.PI) / 180;

  // Convert lat/lng to Cartesian coordinates on unit sphere
  const x1 = Math.cos(lat1) * Math.cos(lng1);
  const y1 = Math.cos(lat1) * Math.sin(lng1);
  const z1 = Math.sin(lat1);

  const x2 = Math.cos(lat2) * Math.cos(lng2);
  const y2 = Math.cos(lat2) * Math.sin(lng2);
  const z2 = Math.sin(lat2);

  // The angle between the two vectors
  let dot = x1 * x2 + y1 * y2 + z1 * z2;
  // Numerical issues can push dot slightly out of [-1,1], so clamp:
  if (dot < -1) dot = -1;
  if (dot > 1) dot = 1;

  const alpha = Math.acos(dot);

  // If alpha=0, the points are the same, no interpolation needed
  if (alpha === 0) {
    return [lat1Deg, lng1Deg];
  }

  const sinAlpha = Math.sin(alpha);
  const ratioA = Math.sin((1 - t) * alpha) / sinAlpha;
  const ratioB = Math.sin(t * alpha) / sinAlpha;

  // Interpolated coordinates in 3D
  const x = ratioA * x1 + ratioB * x2;
  const y = ratioA * y1 + ratioB * y2;
  const z = ratioA * z1 + ratioB * z2;

  // Convert back to spherical
  const interpLat = Math.atan2(z, Math.sqrt(x * x + y * y));
  const interpLng = Math.atan2(y, x);

  // Convert radians back to degrees
  const latDeg = (interpLat * 180) / Math.PI;
  const lngDeg = (interpLng * 180) / Math.PI;

  return [latDeg, lngDeg];
}

/**
 * Linearly interpolate coordinates for the given `time`
 * between `prev` and `next`.
 */
function interpolateCoordinates(prev: Point, next: Point, time: Date): [number, number] {
  // Assert that we're not dealing with an illogical timeline
  console.assert(
    prev.timestamp.getTime() <= next.timestamp.getTime(),
    `prev.timestamp (${prev.timestamp.toISOString()}) must be <= next.timestamp (${next.timestamp.toISOString()})`
  );

  const totalMs = next.timestamp.getTime() - prev.timestamp.getTime();
  const elapsedMs = time.getTime() - prev.timestamp.getTime();

  // Edge case: if both timestamps are identical, just return prev's coords
  if (totalMs === 0) {
    return prev.coordinates;
  }

  // Ensure we aren't interpolating outside the bracket [prev.timestamp, next.timestamp]
  console.assert(
    elapsedMs >= 0 && elapsedMs <= totalMs,
    `Time (${time.toISOString()}) is out of range for interpolation between ${prev.timestamp.toISOString()} and ${next.timestamp.toISOString()}`
  );

  const ratio = elapsedMs / totalMs;

  // Use geodesic interpolation instead of naive linear:
  const [prevLat, prevLng] = prev.coordinates;
  const [nextLat, nextLng] = next.coordinates;

  return interpolateGeodesic(prevLat, prevLng, nextLat, nextLng, ratio);
}

/**
 * Generate a new array of points at 1-minute intervals,
 * from the earliest timestamp to the latest timestamp.
 */
function generateInterpolatedPoints(points: Point[]): Point[] {
  // Assert that we have some points
  console.assert(points.length > 0, "Input array 'points' must not be empty");

  // Sort points by timestamp if not already sorted.
  points.sort((a, b) => a.timestamp.getTime() - b.timestamp.getTime());

  // Re-verify that points are in ascending order
  for (let idx = 0; idx < points.length - 1; idx++) {
    console.assert(
      points[idx].timestamp.getTime() <= points[idx + 1].timestamp.getTime(),
      `Points are not sorted or have duplicates at index ${idx}`
    );
  }

  // If there's only one point, we can't interpolate between two points,
  // so we simply return that single point for every minute in the range.
  if (points.length === 1) {
    // Just replicate the single point's coordinates across whatever range you decide.
    // For example, you might return only that point:
    return [points[0]];
  }

  const startTime = points[0].timestamp;
  const endTime = points[points.length - 1].timestamp;

  // Assert that start <= end
  console.assert(
    startTime.getTime() <= endTime.getTime(),
    "Start time must be <= end time"
  );

  const interpolated: Point[] = [];

  // Pointer for the "previous" point
  let i = 0;

  // 1) Walk minute-by-minute from start to end
  let currentTime = new Date(startTime);
  while (currentTime <= endTime) {
    // 2) If the current time is beyond points[i+1], then move to the next pair
    //    But make sure we don't run off the end of the array.
    while (
      i < points.length - 1 &&
      currentTime.getTime() > points[i + 1].timestamp.getTime()
    ) {
      i++;
    }

    // Edge case: if we're at the very last point, just use its coords
    if (i === points.length - 1) {
      interpolated.push({
        timestamp: new Date(currentTime),
        coordinates: points[i].coordinates,
      });
    } else {
      // Ensure i and i+1 are valid
      console.assert(i + 1 < points.length, "Index out of range for interpolation");

      // 3) Interpolate between points[i] (prev) and points[i+1] (next)
      const prevPoint = points[i];
      const nextPoint = points[i + 1];

      // Make sure currentTime is within or on the edges of [prev, next]
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

    // Advance time by 5 minutes
    currentTime = addMinutes(currentTime, 5);
  }

  return interpolated;
}

type Frame = {
  start: Date;
  end: Date;
  encodedLines: string[];
  center: [number, number];
  zoom: number;
};

/**
 * 1) Weighted lat/lng center in degrees.
 *    (Alternatively, you could use your existing weightedGeodesicCenter
 *     for the center lat & lng.)
 */
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

/**
 * 2) Compute weighted lat/lng standard deviation in kilometers.
 *    RMS => sqrt( E[(Δ km)^2] ).
 *
 *    - For ΔLat: multiply ΔLat(deg) by ~111 km/deg.
 *    - For ΔLng: multiply ΔLng(deg) by ~111 * cos(centerLat) km/deg.
 */
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

  // Standard deviation = sqrt(variance)
  const latStdKm = Math.sqrt(latVariance);
  const lngStdKm = Math.sqrt(lngVariance);

  return { latStdKm, lngStdKm };
}

/**
 * 3) Suppose you want to compute the “max spread” among lat/lng deviations
 *    and pass that to your existing ‘calculateZoomFromDistance()’.
 */
function computeZoomFromLatLngStdDevs(
  weightedPoints: { lat: number; lng: number; weight: number }[]
): number {
  // 1) Find center
  const center = weightedLatLngCenter(weightedPoints);
  // (Or use your existing weightedGeodesicCenter, whichever you prefer.)

  // 2) Find lat/lng standard deviations in km
  const { latStdKm, lngStdKm } = computeWeightedLatLngStdDevKm(weightedPoints, center);

  const aspectRatio = width / height;

  // 3) Use the max
  const maxStdKm = Math.max(latStdKm * aspectRatio, lngStdKm);

  // 4) Pass that to your “calculateZoomFromDistance”
  return calculateZoomFromDistance(maxStdKm);
}

/**
 * Your existing method that takes a distance (km) -> zoom.
 * Very rough mapping, but good enough for illustration.
 */
function calculateZoomFromDistance(distKm: number): number {
  if (distKm <= 0) {
    return 6; // default
  }

  const desiredSpanKm = 4 * distKm; // ±2σ in each direction
  const fraction = 40000 / desiredSpanKm;
  const pxRatio = width / 256; // how many 256-tiles fit horizontally
  let zoom = Math.log2(fraction * pxRatio);

  // clamp & tweak
  zoom = Math.max(3, Math.min(11.8, zoom / 1.4 + 0.4));
  return zoom;
}

/**
 * Weighted geodesic center on the unit sphere.
 * @param weightedPoints An array of objects, each with
 *   { lat: number (deg), lng: number (deg), weight: number }
 * @returns [lat, lng] in degrees
 */
function weightedGeodesicCenter(
  weightedPoints: { lat: number; lng: number; weight: number }[]
): [number, number] {
  let xSum = 0;
  let ySum = 0;
  let zSum = 0;
  let totalWeight = 0;

  for (const { lat, lng, weight } of weightedPoints) {
    // Convert lat/lng to radians
    const latRad = (lat * Math.PI) / 180;
    const lngRad = (lng * Math.PI) / 180;
    // Cartesian on unit sphere
    const cosLat = Math.cos(latRad);
    const x = cosLat * Math.cos(lngRad);
    const y = cosLat * Math.sin(lngRad);
    const z = Math.sin(latRad);
    // Accumulate weighted sums
    xSum += x * weight;
    ySum += y * weight;
    zSum += z * weight;
    totalWeight += weight;
  }

  if (totalWeight === 0) {
    // Edge case: all zero weights
    return [0, 0];
  }

  // Average
  xSum /= totalWeight;
  ySum /= totalWeight;
  zSum /= totalWeight;

  // Convert back to lat/lng
  const hyp = Math.sqrt(xSum * xSum + ySum * ySum);
  const centerLat = Math.atan2(zSum, hyp) * (180 / Math.PI);
  const centerLng = Math.atan2(ySum, xSum) * (180 / Math.PI);

  return [centerLat, centerLng];
}

/**
 * Splits a list of coordinates so that each sub-sequence never crosses ±180.
 * Specifically, if we detect a jump from >+175 to <−175 or vice versa,
 * we “break” the path there, and start a new path from the next point.
 *
 * @param coords Array of [lat, lng]
 * @returns An array of segments, each of which does not cross ±180
 */
function breakAntimeridianSimple(coords: [number, number][]): Array<[number, number][]> {
  if (coords.length <= 1) {
    return [coords];
  }

  const segments: Array<[number, number][]> = [];
  let currentSegment: [number, number][] = [coords[0]];

  for (let i = 1; i < coords.length; i++) {
    const [latPrev, lngPrev] = coords[i - 1];
    const [latCurr, lngCurr] = coords[i];

    // If we jump from >175° to <−175°, or <−175° to >175°, consider it a crossing
    if ((lngPrev > 175 && lngCurr < -175) || (lngPrev < -175 && lngCurr > 175)) {
      // End the current segment at coords[i-1] (already there)
      segments.push(currentSegment);

      // Start a new segment from coords[i]
      currentSegment = [[latCurr, lngCurr]];
    } else {
      // Add the next point to the current segment
      currentSegment.push([latCurr, lngCurr]);
    }
  }

  // Finally push the last segment if it has any coords
  if (currentSegment.length > 0) {
    segments.push(currentSegment);
  }

  return segments;
}

function generateFrames(
  points: Point[],
  animationStart: Date,
  animationEnd: Date,
  stepInHours: number = 1,
  pathWindowInHours: number = 24,
  centerWindowInHours: number = 120
) {
  // Prefix sums for faster path-encoding (unchanged)
  const latPrefix = new Array<number>(points.length + 1).fill(0);
  const lngPrefix = new Array<number>(points.length + 1).fill(0);
  for (let i = 0; i < points.length; i++) {
    latPrefix[i + 1] = latPrefix[i] + points[i].coordinates[0];
    lngPrefix[i + 1] = lngPrefix[i] + points[i].coordinates[1];
  }

  const frames: Array<Frame> = [];

  // For slicing the path (24h window)
  let pathLeft = 0;
  let pathRight = 0;

  // For slicing the center (120h window)
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

    // Path window: [mid - 12h, mid + 12h]
    const pathStartMs = frameMidMs - pathWindowSide;
    const pathEndMs = frameMidMs + pathWindowSide;

    // Move pathLeft to first point >= pathStartMs
    while (
      pathLeft < points.length &&
      points[pathLeft].timestamp.getTime() < pathStartMs
    ) {
      pathLeft++;
    }
    // Move pathRight to first point > pathEndMs
    while (
      pathRight < points.length &&
      points[pathRight].timestamp.getTime() <= pathEndMs
    ) {
      pathRight++;
    }

    const pathPoints = points.slice(pathLeft, pathRight);
    const pathCoords = pathPoints.map((p) => p.coordinates);

    const segments = breakAntimeridianSimple(pathCoords);

    // Then encode each segment independently:
    const encodedLines = segments.map((segment) => encode(segment));

    // ---- Weighted center calculation ----
    // For the center (120h window): [mid - 60h, mid + 60h]
    const centerStartMs = frameMidMs - centerWindowSide;
    const centerEndMs = frameMidMs + centerWindowSide;

    // Move centerLeft to first point >= centerStartMs
    while (
      centerLeft < points.length &&
      points[centerLeft].timestamp.getTime() < centerStartMs
    ) {
      centerLeft++;
    }
    // Move centerRight to first point > centerEndMs
    while (
      centerRight < points.length &&
      points[centerRight].timestamp.getTime() <= centerEndMs
    ) {
      centerRight++;
    }

    // Triangular weighting function around frameMidMs
    const halfWindow = centerWindowSide;

    function gaussianWeight(timestampMs: number, centerMs: number) {
      const dist = Math.abs(timestampMs - centerMs);
      // Pick some sigma ~ (halfWindow / 2) or whichever you prefer
      const sigma = halfWindow / 2;
      return Math.exp(-0.5 * (dist / sigma) ** 2);
    }

    // Weighted sums in 3D
    const weightedPoints = [];

    for (let idx = centerLeft; idx < centerRight; idx++) {
      const { timestamp, coordinates } = points[idx];
      const [lat, lng] = coordinates;
      const w = gaussianWeight(timestamp.getTime(), frameMidMs);
      if (w > 0) {
        // Store for geodesic center
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

/**
 * Smooth the `zoom` values of each frame by applying a Gaussian kernel
 * that looks 2 frames before and 2 frames after.
 *
 * @param frames The array of frames, each of which has a `zoom` property to be smoothed.
 * @param radius The number of frames on each side of the current frame to include.
 * @param sigma The standard deviation of the Gaussian kernel.
 */
function smoothZooms(frames: Frame[], radius = 2, sigma = 1.0) {
  // 1) Generate the kernel weights
  const kernelSize = radius * 2 + 1; // e.g. radius=2 => kernelSize=5
  const kernel = new Array<number>(kernelSize);

  let sum = 0;
  for (let i = -radius; i <= radius; i++) {
    const x = (i * i) / (2 * sigma * sigma);
    // Basic Gaussian kernel: e^{-x}
    const w = Math.exp(-x);
    kernel[i + radius] = w;
    sum += w;
  }

  // Normalize the kernel so that all weights sum to 1
  for (let i = 0; i < kernelSize; i++) {
    kernel[i] /= sum;
  }

  // 2) Create an array to hold the smoothed zoom values
  const smoothedZooms = new Array<number>(frames.length);

  // 3) For each frame, apply the kernel
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

  // 4) Store the smoothed zoom values back into frames
  for (let i = 0; i < frames.length; i++) {
    frames[i].zoom = smoothedZooms[i];
  }
}

const interpolatedPoints = generateInterpolatedPoints(points);

const frames = generateFrames(interpolatedPoints, animationStart, animationEnd);

console.log("Generated", frames.length, "frames-test");

// Smooth the zoom values

smoothZooms(frames, 168, 15.0);

console.log("Smoothed zoom values");

// build url for each frame
const accessToken =
  "pk.eyJ1Ijoiam9zaGNoYW5nMDQiLCJhIjoiY2p3c2c5NDdsMDEyOTQwcXhyc245eTYwcSJ9.CKHrw2ddbsyG9bdOlr0TvQ";
const stroke = "ff0000";
const weight = 3;

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

// fetch & save
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
