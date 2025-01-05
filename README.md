# animate-timeline

This project generates an animated video of your location history from a Google Maps Timeline JSON file.

## Usage

Clone this repository:

```bash
git clone https://github.com/joshdchang/animate-timeline.git
cd animate-timeline
```

You will need bun to run this project. To install bun:

```bash
curl -fsSL https://bun.sh/install | bash
```

You will also need ffmpeg to generate the video. To install ffmpeg on macOS:

```bash
brew install ffmpeg
```

Then, to install dependencies:

```bash
bun install
```

To run:

```bash
bun run index.ts <location-history.json>
```

Where `<location-history.json>` is the path to the location history JSON file exported from your Google Maps Timeline on your phone. You can find this by clicking on your profile picture in Google Maps, then "Your Timeline", then the three dots in the top right, then "Location & privacy settings", then scroll down and click "Export Timeline data" and send it to your computer.
