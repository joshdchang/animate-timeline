import sharp from "sharp";
import fs from "fs/promises";
import { animationStart } from "./config";

await fs.mkdir("text-frames", { recursive: true });

const frames = await fs.readdir("frames");

for (const frame of frames) {
  const framePath = `frames/${frame}`;
  const frameImage = sharp(framePath);

  const frameNumber = parseInt(frame.split("-")[1].split(".")[0]);

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
          dpi: 300, // dots per inch
          rgba: true,
        },
      },
      blend: "darken",
      gravity: "northwest",
      left: 40,
      top: 30,
    },
  ]);

  frameWithText.toFile(`text-frames/${frame}`).then(() => {
    console.log(`Frame ${frame} processed`);
  });
}
