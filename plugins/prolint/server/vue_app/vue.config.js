const path = require("path");

module.exports = {
  outputDir: path.resolve(__dirname, "../dist/"),
  assetsDir: "./assets",
  transpileDependencies: ["vuetify"],
};
